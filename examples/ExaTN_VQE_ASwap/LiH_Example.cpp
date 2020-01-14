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
    // LiH (at bond distance) Hamiltonian from https://arxiv.org/abs/1704.05018 (Supplementary information)
    auto LiH = xacc::quantum::getObservable(
      "pauli", std::string(R"(
            -0.096022 Z0 + 
            -0.206128 Z0Z1 + 
            0.364746 Z1 + 
            0.096022 Z2 + 
            -0.206128 Z2Z3 + 
            -0.364746 Z3 + 
            -0.145438 Z0Z2 + 
            0.056040 Z0Z2Z3 + 
            0.110811 Z0Z3 + 
            -0.056040 Z0Z1Z2 + 
            0.080334 Z0Z1Z2Z3 + 
            0.063673 Z0Z1Z3 + 
            0.110811 Z1Z2 + 
            -0.063673 Z1Z2Z3 + 
            -0.095216 Z1Z3 + 

            -0.012585 X0Z1 + 
            0.012585 X0 + 
            0.012585 X2Z3 + 
            0.012585 X2 + 
            -0.002667 X0Z1X2Z3 + 
            -0.002667 X0Z1X2 + 
            0.002667 X0X2Z3 + 
            0.002667 X0X2 + 
            0.007265 X0Z1Z3 + 
            -0.007265 X0Z3 + 
            0.007265 Z1X2Z3 + 
            0.007265 Z1X2 + 


            -0.029640 X0X1 + 
            0.002792 X1 + 
            -0.029640 X2X3 + 
            0.002792 X3 + 
            -0.008195 X0X2X3 + 
            -0.001271 X0X3 + 
            -0.008195 X0X1X2 + 
            0.028926 X0X1X2X3 + 
            0.007499 X0X1X3 + 
            -0.001271 X1X2 + 
            0.007499 X1X2X3 + 
            0.009327 X1X3 + 

            0.029640 Y0Y1 + 
            0.029640 Y2Y3 + 
            0.028926 Y0Y1Y2Y3 + 

            0.002792 Z0X1 + 
            -0.002792 Z2X3 + 
            -0.016781 Z0Z2X3 + 
            0.016781 Z0X3 + 
            -0.016781 Z0X1Z2 + 
            -0.016781 X1Z2 + 
            -0.009327 Z0X1Z2X3 + 
            0.009327 Z0X1X3 + 
            -0.009327 X1Z2X3 + 

            -0.011962 Z0X2Z3 + 
            -0.011962 Z0X2 + 
            0.000247 Z0Z1X2Z3 + 
            0.000247 Z0Z1X2 + 

            0.039155 Z0X2X3 + 
            -0.002895 Z0Z1X2X3 + 
            -0.009769 Z0Z1X3 + 
            -0.024280 Z1X2X3 + 
            -0.008025 Z1X3 + 

            -0.039155 Z0Y2Y3 + 
            0.002895 Z0Z1Y2Y3 + 
            0.024280 Z1Y2Y3 + 

            -0.011962 X0Z1Z2 + 
            0.011962 X0Z2 + 
            -0.000247 X0Z1Z2Z3 + 
            0.000247 X0Z2Z3 + 

            0.008195 X0Z1X2X3 + 
            0.001271 X0Z1X3 + 

            -0.008195 X0Z1Y2Y3 + 
            0.008195 X0Y2Y3 + 

            -0.001271 X0Z1Z2X3 + 
            0.001271 X0Z2X3 + 
            0.008025 Z1Z2X3 + 

            -0.039155 X0X1Z2 + 
            -0.002895 X0X1Z2Z3 + 
            0.024280 X0X1Z3 + 
            -0.009769 X1Z2Z3 + 
            0.008025 X1Z3 + 

            0.039155 Y0Y1Z2 + 
            0.002895 Y0Y1Z2Z3 + 
            -0.024280 Y0Y1Z3 + 

            -0.008195 X0X1X2Z3 + 
            -0.001271 X1X2Z3 + 

            0.008195 Y0Y1X2Z3 + 
            0.008195 Y0Y1X2 + 

            -0.028926 X0X1Y2Y3 + 
            -0.007499 X1Y2Y3 + 

            -0.028926 Y0Y1X2X3 + 
            -0.007499 Y0Y1X3 + 

            -0.007499 X0X1Z2X3 + 

            0.007499 Y0Y1Z2X3 + 

            0.009769 Z0Z1Z2X3 + 

            -0.001271 Z0X1X2Z3 + 
            -0.001271 Z0X1X2 + 
            0.008025 Z0X1Z3 + 

            0.007499 Z0X1X2X3 + 

            -0.007499 Z0X1Y2Y3 + 

            -0.009769 Z0X1Z2Z3 
    )"));

    const auto pauliObservable = std::static_pointer_cast<xacc::quantum::PauliOperator>(LiH);    
    const auto tnqvmAccelerator = std::static_pointer_cast<tnqvm::TNQVM>(qpu);

    // ASWAP with 4 qubits (orbitals) and 2 particles (electron)
    xacc::qasm(R"(
    .compiler xasm
    .circuit LiH_ansatz
    .parameters t0, t1, t2, t3, t4, t5
    .qbit q
    ASWAP(q, t0, t1, t2, t3, t4, t5, {{"nbQubits", 4}, {"nbParticles", 2}});
    )");

    auto ansatz = xacc::getCompiled("LiH_ansatz");

    // Debug:
    std::cout << "Circuit: \n" << ansatz->toString();

    auto buffer = xacc::qalloc(4);
    
    size_t trialCount = 0;

    xacc::OptFunction optFunc([&](const std::vector<double>& x, std::vector<double>& dx) {
        const auto start = std::chrono::steady_clock::now();
        
        auto evaled = ansatz->operator()(x);
        const double expVal = tnqvmAccelerator->getExpectationValue(buffer, evaled, *pauliObservable);
        const auto end = std::chrono::steady_clock::now();

        // Logging
        {
            std::cout << "Trial#" << trialCount << ": Elapsed time: " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << " ms \n";
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
                    << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                    << " ms \n";
    }
    
    buffer->addExtraInfo("opt-val", ExtraInfo(result.first));
    buffer->addExtraInfo("opt-params", ExtraInfo(result.second));
     
    std::cout << "Energy: " << buffer->getInformation("opt-val").as<double>() << "\n";
    
    // Finalize the XACC Framework
	xacc::Finalize();

	return 0;
}
