/***********************************************************************************
 * Copyright (c) 2020, UT-Battelle
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the xacc nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * Contributors:
 *   Initial implementation - Thien Nguyen
 * 
**********************************************************************************/
#include "xacc.hpp"
#include "TNQVM.hpp"
#include <chrono>
#include <unistd.h>
#include "xacc_observable.hpp"
#include "PauliOperator.hpp"

int main (int argc, char** argv) {
    using namespace tnqvm;

	// Initialize the XACC Framework
	xacc::Initialize(argc, argv);
       
	auto qpu = xacc::getAccelerator("tnqvm", { std::make_pair("tnqvm-visitor", "exatn") });
	auto optimizer = xacc::getOptimizer("nlopt");
    
    const auto tnqvmAccelerator = std::static_pointer_cast<tnqvm::TNQVM>(qpu);   
    auto H_N_3 = xacc::quantum::getObservable(
      "pauli",
      std::string("5.907 - 2.1433 X0X1 - 2.1433 Y0Y1 + .21829 Z0 - 6.125 Z1 + "
                  "9.625 - 9.625 Z2 - 3.91 X1 X2 - 3.91 Y1 Y2"));
    const auto pauliObservable = std::static_pointer_cast<xacc::quantum::PauliOperator>(H_N_3);
    
    // ASWAP with 3 qubits (orbitals) and 1 particle (electron)
    xacc::qasm(R"(
    .compiler xasm
    .circuit deuteron_ansatz_h3
    .parameters t0, t1, t2
    .qbit q
    ASWAP(q, t0, t1, t2, {{"nbQubits", 3}, {"nbParticles", 1}});
    )");

    auto ansatz = xacc::getCompiled("deuteron_ansatz_h3");
    auto buffer = xacc::qalloc(3);
    size_t trialCount = 0;

    xacc::OptFunction optFunc([&](const std::vector<double>& x, std::vector<double>& dx) {
        const auto start = std::chrono::steady_clock::now();
        auto evaled = ansatz->operator()(x);
        const double expVal = tnqvmAccelerator->getExpectationValue(buffer, evaled, *pauliObservable);
        const auto end = std::chrono::steady_clock::now();

        // Logging
        {
            std::cout << "Trial#" << trialCount << ": Elapsed time in seconds : " 
                        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                        << " secs \n";

            std::cout << "E(" << x[0];
            for (int i = 1; i < x.size(); ++i)
            {
                std::cout << "," << x[i];
            }          
            std::cout << ") = " << expVal << "\n";
        }
                
        trialCount++;
        return expVal;
    },
    ansatz->nVariables()
    );

    const auto start = std::chrono::steady_clock::now();
    auto result = optimizer->optimize(optFunc);
    const auto end = std::chrono::steady_clock::now();


    {
        std::cout << "Complete: number of trials = " << trialCount << ", total elapsed time = "
                    << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                    << " secs \n";
    }
    
    buffer->addExtraInfo("opt-val", ExtraInfo(result.first));
    buffer->addExtraInfo("opt-params", ExtraInfo(result.second));
     
    // Expected result: -2.04482
    std::cout << "Energy: " << buffer->getInformation("opt-val").as<double>() << "\n"; 
  
	// Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}
