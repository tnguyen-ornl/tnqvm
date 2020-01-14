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
    // BeH2 (at bond distance) Hamiltonian from https://arxiv.org/abs/1704.05018 (Supplementary information)
    auto BeH2 = xacc::quantum::getObservable(
      "pauli", std::string(R"(
            -0.143021 Z0 + 
            0.104962 Z0Z1 + 
            0.038195 Z1Z2 + 
            -0.325651 Z2 + 
            -0.143021 Z3 + 
            0.104962 Z3Z4 + 
            0.038195 Z4Z5 + 
            -0.325651 Z5 + 
            0.172191 Z1 + 
            0.174763 Z0Z1Z2 + 
            0.136055 Z0Z2 + 
            0.116134 Z0Z3 + 
            0.094064 Z0Z3Z4 + 
            0.099152 Z0Z4Z5 + 
            0.123367 Z0Z5 + 
            0.094064 Z0Z1Z3 + 
            0.098003 Z0Z1Z3Z4 + 
            0.102525 Z0Z1Z4Z5 + 
            0.097795 Z0Z1Z5 + 
            0.099152 Z1Z2Z3 + 
            0.102525 Z1Z2Z3Z4 + 
            0.112045 Z1Z2Z4Z5 + 
            0.105708 Z1Z2Z5 + 
            0.123367 Z2Z3 + 
            0.097795 Z2Z3Z4 + 
            0.105708 Z2Z4Z5 + 
            0.133557 Z2Z5 + 
            0.172191 Z4 + 
            0.174763 Z3Z4Z5 + 
            0.136055 Z3Z5 + 
            0.059110 X0Z1 + 
            -0.059110 X0 + 
            0.161019 Z1X2 + 
            -0.161019 X2 + 
            0.059110 X3Z4 + 
            -0.059110 X3 + 
            0.161019 Z4X5 + 
            -0.161019 X5 + 
            -0.038098 X0X2 + 
            -0.003300 X0Z1X2 + 
            0.013745 X0Z1X3Z4 + 
            -0.013745 X0Z1X3 + 
            -0.013745 X0X3Z4 + 
            0.013745 X0X3 + 
            0.011986 X0Z1Z4X5 + 
            -0.011986 X0Z1X5 + 
            -0.011986 X0Z4X5 + 
            0.011986 X0X5 + 
            0.011986 Z1X2X3Z4 + 
            -0.011986 Z1X2X3 + 
            -0.011986 X2X3Z4 + 
            0.011986 X2X3 + 
            0.013836 Z1X2Z4X5 + 
            -0.013836 Z1X2X5 + 
            -0.013836 X2Z4X5 + 
            0.013836 X2X5 + 
            -0.038098 X3X5 + 
            -0.003300 X3Z4X5 + 
            -0.002246 Z0Z1X2 + 
            0.002246 Z0X2 + 
            0.014815 Z0X3Z4 + 
            -0.014815 Z0X3 + 
            0.009922 Z0Z4X5 + 
            -0.009922 Z0X5 + 
            -0.002038 Z0Z1X3Z4 + 
            0.002038 Z0Z1X3 + 
            -0.007016 Z0Z1Z4X5 + 
            0.007016 Z0Z1X5 + 
            -0.006154 X0Z2 + 
            0.006154 X0Z1Z2 + 
            0.014815 X0Z1Z3 + 
            -0.014815 X0Z3 + 
            -0.002038 X0Z1Z3Z4 + 
            0.002038 X0Z3Z4 + 
            0.001124 X0Z1Z4Z5 + 
            -0.001124 X0Z4Z5 + 
            0.017678 X0Z1Z5 + 
            -0.017678 X0Z5 + 
            -0.041398 Y0Y2 + 
            0.011583 Y0Y1X3X4Z5 + 
            -0.011094 Y0Y1X4 + 
            0.010336 Y1Y2X3X4Z5 + 
            -0.005725 Y1Y2X4 + 
            -0.006154 X3Z5 + 
            0.011583 X0X1Z2X3X4Z5 + 
            -0.011094 X0X1Z2X4 + 
            -0.011094 X1X3X4Z5 + 
            0.026631 X1X4 + 
            -0.017678 Z2X3 + 
            0.011583 X0X1Z2Y3Y4 + 
            0.010336 X0X1Z2Y4Y5 + 
            -0.011094 X1Y3Y4 + 
            -0.005725 X1Y4Y5 + 
            -0.041398 Y3Y5 + 
            0.011583 Y0Y1Y3Y4 + 
            0.010336 Y0Y1Y4Y5 + 
            0.010336 Y1Y2Y3Y4 + 
            0.010600 Y1Y2Y4Y5 + 
            0.024909 X0X1Z2X3X4X5 + 
            -0.031035 X1X3X4X5 + 
            -0.010064 Z2X5 + 
            0.024909 X0X1Z2Y3X4Y5 + 
            -0.031035 X1Y3X4Y5 + 
            0.024909 Y0Y1X3X4X5 + 
            0.021494 Y1Y2X3X4X5 + 
            0.024909 Y0Y1Y3X4Y5 + 
            0.021494 Y1Y2Y3X4Y5 + 
            0.011094 X0X1Z2Z3X4Z5 + 
            -0.026631 X1Z3X4Z5 + 
            0.011094 Y0Y1Z3X4Z5 + 
            0.005725 Y1Y2Z3X4Z5 + 
            0.010336 X0X1Z2Z3X4X5 + 
            -0.005725 X1Z3X4X5 + 
            0.002246 Z3X5 + 
            0.010336 Y0Y1Z3X4X5 + 
            0.010600 Y1Y2Z3X4X5 + 
            0.024909 X0X1X2X3X4Z5 + 
            -0.031035 X0X1X2X4 + 
            -0.010064 X2Z5 + 
            0.024909 X0X1X2Y3Y4 + 
            0.021494 X0X1X2Y4Y5 + 
            0.024909 Y0X1Y2X3X4Z5 + 
            -0.031035 Y0X1Y2X4 + 
            0.024909 Y0X1Y2Y3Y4 + 
            0.021494 Y0X1Y2Y4Y5 + 
            0.063207 X0X1X2X3X4X5 + 
            0.063207 X0X1X2Y3X4Y5 + 
            0.063207 Y0X1Y2X3X4X5 + 
            0.063207 Y0X1Y2Y3X4Y5 + 
            0.031035 X0X1X2Z3X4Z5 + 
            -0.009922 X2Z3 + 
            0.031035 Y0X1Y2Z3X4Z5 + 
            0.021494 X0X1X2Z3X4X5 + 
            0.021494 Y0X1Y2Z3X4X5 + 
            0.011094 Z0X1Z2X3X4Z5 + 
            -0.026631 Z0X1Z2X4 + 
            0.011094 Z0X1Z2Y3Y4 + 
            0.005725 Z0X1Z2Y4Y5 + 
            0.031035 Z0X1Z2X3X4X5 + 
            0.031035 Z0X1Z2Y3X4Y5 + 
            0.026631 Z0X1Z2Z3X4Z5 + 
            0.005725 Z0X1Z2Z3X4X5 + 
            0.010336 Z0X1X2X3X4Z5 + 
            -0.005725 Z0X1X2X4 + 
            0.010336 Z0X1X2Y3Y4 + 
            0.010600 Z0X1X2Y4Y5 + 
            0.021494 Z0X1X2X3X4X5 + 
            0.021494 Z0X1X2Y3X4Y5 + 
            0.005725 Z0X1X2Z3X4Z5 + 
            0.010600 Z0X1X2Z3X4X5 + 
            0.001124 Z1Z2X3Z4 + 
            -0.001124 Z1Z2X3 + 
            -0.007952 Z1Z2Z4X5 + 
            0.007952 Z1Z2X5 + 
            0.017678 Z2X3Z4 + 
            0.010064 Z2Z4X5 + 
            0.009922 Z1X2Z3 + 
            -0.007016 Z1X2Z3Z4 + 
            0.007016 X2Z3Z4 + 
            -0.007952 Z1X2Z4Z5 + 
            0.007952 X2Z4Z5 + 
            0.010064 Z1X2Z5 + 
            -0.002246 Z3Z4X5 + 
            0.006154 X3Z4Z5
    )"));

    const auto pauliObservable = std::static_pointer_cast<xacc::quantum::PauliOperator>(BeH2);    
    const auto tnqvmAccelerator = std::static_pointer_cast<tnqvm::TNQVM>(qpu);

    // ASWAP with 6 qubits (orbitals) and 4 particles (electron)
    xacc::qasm(R"(
    .compiler xasm
    .circuit BeH2_ansatz
    .parameters t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14
    .qbit q
    ASWAP(q, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, {{"nbQubits", 6}, {"nbParticles", 4}});
    )");

    auto ansatz = xacc::getCompiled("BeH2_ansatz");

    // Debug:
    std::cout << "Circuit: \n" << ansatz->toString();

    auto buffer = xacc::qalloc(6);
    
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
