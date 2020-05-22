#include "xacc.hpp"
#include "xacc_service.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include <chrono>

int main(int argc, char **argv) 
{
    xacc::Initialize();

    auto tmp = xacc::getService<xacc::Instruction>("rcs");
    auto randomCirc = std::dynamic_pointer_cast<xacc::CompositeInstruction>(tmp);
    const int NB_QUBITS = 21;
    randomCirc->expand({
        std::make_pair("nq", NB_QUBITS), 
        std::make_pair("nlayers", 8),
        std::make_pair("parametric-gates", false)
    });

    // Number of shots must be large enough to observe the distribution.
    const int NB_SHOTS = 1 << NB_QUBITS;
    // We'll save the output to file for processing.
    std::ofstream out("output.txt");
    out << "Number of gates = " << randomCirc->nInstructions() << "\n";
    out << "Circuit:\n" << randomCirc->toString() << "\n";

    auto accelerator = xacc::getAccelerator("tnqvm", {
        std::make_pair("tnqvm-visitor", "exatn"), 
        std::make_pair("shots", NB_SHOTS),
    });

    auto qreg = xacc::qalloc(NB_QUBITS);
    const auto start = std::chrono::system_clock::now();
    accelerator->execute(qreg, randomCirc);
    const auto end = std::chrono::system_clock::now();
    out << "Elapsed time in milliseconds : " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
    << " ms\n";
    qreg->print();
    // Save to file.
    out << qreg->toString();
    xacc::Finalize();
    
    out.close();
    return 0;
} 
