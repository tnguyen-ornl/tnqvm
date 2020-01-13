add_executable(DeuteuronH2Example DeuteuronExampleH2.cpp)
target_link_libraries(DeuteuronH2Example PRIVATE xacc::xacc tnqvm)

add_executable(DeuteuronH3Example DeuteuronExampleH3.cpp)
target_link_libraries(DeuteuronH3Example PRIVATE xacc::xacc tnqvm)