#include <memory>
#include <gtest/gtest.h>
#include "xacc.hpp"
#include "xacc_service.hpp"

using namespace xacc;

namespace {
const std::string srcStr = R"(
__qpu__ void sycamoreCirc(qbit q) {
// Begin hz_1_2
Rz(q[0], -0.78539816339);
Rx(q[0], 1.57079632679);
Rz(q[0], 0.78539816339);
// End hz_1_2
Rx(q[1], 1.57079632679);
Rx(q[2], 1.57079632679);
// Begin hz_1_2
Rz(q[3], -0.78539816339);
Rx(q[3], 1.57079632679);
Rz(q[3], 0.78539816339);
// End hz_1_2
Ry(q[4], 1.57079632679);
// Begin hz_1_2
Rz(q[5], -0.78539816339);
Rx(q[5], 1.57079632679);
Rz(q[5], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[6], -0.78539816339);
Rx(q[6], 1.57079632679);
Rz(q[6], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[7], -0.78539816339);
Rx(q[7], 1.57079632679);
Rz(q[7], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[8], -0.78539816339);
Rx(q[8], 1.57079632679);
Rz(q[8], 0.78539816339);
// End hz_1_2
Rx(q[9], 1.57079632679);
Ry(q[10], 1.57079632679);
// Begin hz_1_2
Rz(q[11], -0.78539816339);
Rx(q[11], 1.57079632679);
Rz(q[11], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[12], -0.78539816339);
Rx(q[12], 1.57079632679);
Rz(q[12], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[13], -0.78539816339);
Rx(q[13], 1.57079632679);
Rz(q[13], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[14], -0.78539816339);
Rx(q[14], 1.57079632679);
Rz(q[14], 0.78539816339);
// End hz_1_2
Ry(q[15], 1.57079632679);
Rx(q[16], 1.57079632679);
Rx(q[17], 1.57079632679);
// Begin hz_1_2
Rz(q[18], -0.78539816339);
Rx(q[18], 1.57079632679);
Rz(q[18], 0.78539816339);
// End hz_1_2
Rx(q[19], 1.57079632679);
Rx(q[20], 1.57079632679);
// Begin hz_1_2
Rz(q[21], -0.78539816339);
Rx(q[21], 1.57079632679);
Rz(q[21], 0.78539816339);
// End hz_1_2
Ry(q[22], 1.57079632679);
Rx(q[23], 1.57079632679);
Rx(q[24], 1.57079632679);
Ry(q[25], 1.57079632679);
// Begin hz_1_2
Rz(q[26], -0.78539816339);
Rx(q[26], 1.57079632679);
Rz(q[26], 0.78539816339);
// End hz_1_2
Rx(q[27], 1.57079632679);
// Begin hz_1_2
Rz(q[28], -0.78539816339);
Rx(q[28], 1.57079632679);
Rz(q[28], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[29], -0.78539816339);
Rx(q[29], 1.57079632679);
Rz(q[29], 0.78539816339);
// End hz_1_2
Rx(q[30], 1.57079632679);
Rx(q[31], 1.57079632679);
Rx(q[32], 1.57079632679);
Ry(q[33], 1.57079632679);
// Begin hz_1_2
Rz(q[34], -0.78539816339);
Rx(q[34], 1.57079632679);
Rz(q[34], 0.78539816339);
// End hz_1_2
Ry(q[35], 1.57079632679);
// Begin hz_1_2
Rz(q[36], -0.78539816339);
Rx(q[36], 1.57079632679);
Rz(q[36], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[37], -0.78539816339);
Rx(q[37], 1.57079632679);
Rz(q[37], 0.78539816339);
// End hz_1_2
// Begin hz_1_2
Rz(q[38], -0.78539816339);
Rx(q[38], 1.57079632679);
Rz(q[38], 0.78539816339);
// End hz_1_2
Ry(q[39], 1.57079632679);
// Begin hz_1_2
Rz(q[40], -0.78539816339);
Rx(q[40], 1.57079632679);
Rz(q[40], 0.78539816339);
// End hz_1_2
Rx(q[41], 1.57079632679);
Rx(q[42], 1.57079632679);
Rx(q[43], 1.57079632679);
// Begin hz_1_2
Rz(q[44], -0.78539816339);
Rx(q[44], 1.57079632679);
Rz(q[44], 0.78539816339);
// End hz_1_2
Ry(q[45], 1.57079632679);
Rx(q[46], 1.57079632679);
// Begin hz_1_2
Rz(q[47], -0.78539816339);
Rx(q[47], 1.57079632679);
Rz(q[47], 0.78539816339);
// End hz_1_2
Ry(q[48], 1.57079632679);
Ry(q[49], 1.57079632679);
Ry(q[50], 1.57079632679);
// Begin hz_1_2
Rz(q[51], -0.78539816339);
Rx(q[51], 1.57079632679);
Rz(q[51], 0.78539816339);
// End hz_1_2
Rx(q[52], 1.57079632679);
Rz(q[1], 2.432656295030221);
Rz(q[4], -2.2258827283782336);
Rz(q[3], -2.7293249642089283);
Rz(q[7], 1.2106965020980647);
Rz(q[5], -1.1065190593715701);
Rz(q[9], 1.7892230375390792);
Rz(q[6], 0.21199582799561975);
Rz(q[13], -0.1112833809595124);
Rz(q[8], 2.893794754566797);
Rz(q[15], -2.954998228479576);
Rz(q[10], 1.28422274116456);
Rz(q[17], -1.1227354153033684);
Rz(q[12], 1.34658993786433);
Rz(q[21], -1.7818294429115995);
Rz(q[14], 2.1872907310718426);
Rz(q[23], -1.961401963214366);
Rz(q[16], 1.5928715721088071);
Rz(q[25], -1.5401880072084049);
Rz(q[18], -2.591169533039924);
Rz(q[27], 2.612251899130507);
Rz(q[20], 2.4047981541022003);
Rz(q[30], -2.394731656273015);
Rz(q[22], -2.3659320009488227);
Rz(q[32], 2.2191901639170335);
Rz(q[24], -2.435037316131151);
Rz(q[34], 3.0221985839801064);
Rz(q[26], -2.6108064207399715);
Rz(q[36], 2.560219630921921);
Rz(q[29], 1.7858344832951236);
Rz(q[37], -1.7147844803366754);
Rz(q[31], 1.4590820918159937);
Rz(q[39], -2.7148644518627782);
Rz(q[33], 1.056801934900941);
Rz(q[41], -1.2947483188429902);
Rz(q[35], 2.768299098244506);
Rz(q[43], -2.576048471247718);
Rz(q[38], -1.6256583901852042);
Rz(q[44], 1.6003028045332237);
Rz(q[40], -0.8433460214813505);
Rz(q[46], 0.8404319695058231);
Rz(q[42], 2.4767110913273673);
Rz(q[48], 2.9232217616730027);
Rz(q[45], 1.7188939810777453);
Rz(q[49], -1.8507490475427473);
Rz(q[47], 2.3200371653685714);
Rz(q[51], -2.3409287132563574);
Rz(q[50], -0.5750038022731371);
Rz(q[52], -0.5519833389562016);
fSim(q[1], q[4], 1.5157741664070026, 0.5567125777724111);
fSim(q[3], q[7], 1.5177580142210796, 0.49481085782254924);
fSim(q[5], q[9], 1.603673862122088, 0.47689957001761957);
fSim(q[6], q[13], 1.517732708964379, 0.5058312223892918);
fSim(q[8], q[15], 1.52531844440771, 0.46557175536522444);
fSim(q[10], q[17], 1.6141004604575337, 0.4943434406753406);
fSim(q[12], q[21], 1.5476810407276227, 0.44290174465705406);
fSim(q[14], q[23], 1.5237261387831182, 0.46966161228464703);
fSim(q[16], q[25], 1.52858449420213, 0.5736654641907271);
fSim(q[18], q[27], 1.5483159975150524, 0.4961408893973949);
fSim(q[20], q[30], 1.6377079485606105, 0.6888985951517979);
fSim(q[22], q[32], 1.5299499142361528, 0.4825884757470581);
fSim(q[24], q[34], 1.5280421758408222, 0.5109767145463228);
fSim(q[26], q[36], 1.512078286877267, 0.48151528098618757);
fSim(q[29], q[37], 1.5071938854286824, 0.5089276063739601);
fSim(q[31], q[39], 1.5460100224552222, 0.5302403303961926);
fSim(q[33], q[41], 1.516662594039817, 0.45171597904279737);
fSim(q[35], q[43], 1.4597689731865275, 0.42149859585364335);
fSim(q[38], q[44], 1.5356494456905732, 0.47076284376184807);
fSim(q[40], q[46], 1.5179778495709582, 0.5221350266177678);
fSim(q[42], q[48], 1.4969321270214224, 0.4326117171327447);
fSim(q[45], q[49], 1.51149872016387, 0.4914319343688027);
fSim(q[47], q[51], 1.4908807480931237, 0.48862437201319);
fSim(q[50], q[52], 1.6162569997269376, 0.5014289362839901);
}
)";
}

// Test Sycamore (RCS) circuits
// Using full tensor network contraction:
TEST(SycamoreCircTester, checkExatn) 
{
    // The bitstring to calculate amplitude
    // Example: bitstring = 000000000...00
    const std::vector<int> BIT_STRING(53, 0);

    // ExaTN visitor: 
    // Note: 
    auto qpu = xacc::getAccelerator("tnqvm", {
        std::make_pair("tnqvm-visitor", "exatn"),
        std::make_pair("bitstring", BIT_STRING)
    });

    // Allocate a register of 53 qubits
    auto qubitReg = xacc::qalloc(53);
    // Compile the XASM program
    auto xasmCompiler = xacc::getCompiler("xasm");
    auto ir = xasmCompiler->compile(srcStr, qpu);
    auto program = ir->getComposites()[0];
    qpu->execute(qubitReg, program);
    // qubitReg->print();
    const double realAmpl = (*qubitReg)["amplitude-real"].as<double>();
    const double imagAmpl = (*qubitReg)["amplitude-imag"].as<double>();
    EXPECT_NEAR(realAmpl, 2.31963e-09, 1e-10);
    EXPECT_NEAR(imagAmpl, -1.02782e-08, 1e-10);
}

int main(int argc, char **argv) 
{
    xacc::Initialize();
    ::testing::InitGoogleTest(&argc, argv);
    auto ret = RUN_ALL_TESTS();
    xacc::Finalize();
    return ret;
} 
