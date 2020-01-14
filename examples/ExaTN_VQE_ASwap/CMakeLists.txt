add_executable(DeuteuronH2Example DeuteuronExampleH2.cpp)
target_link_libraries(DeuteuronH2Example PRIVATE xacc::xacc tnqvm)

add_executable(DeuteuronH3Example DeuteuronExampleH3.cpp)
target_link_libraries(DeuteuronH3Example PRIVATE xacc::xacc tnqvm)

add_executable(LiH_Example LiH_Example.cpp)
target_link_libraries(LiH_Example PRIVATE xacc::xacc tnqvm)

add_executable(BeH2_Example BeH2_Example.cpp)
target_link_libraries(BeH2_Example PRIVATE xacc::xacc tnqvm)