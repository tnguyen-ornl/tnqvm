include_directories(${XACC_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR})

add_executable(check_bitstring_distribution check_bitstring_distribution.cpp)
target_link_libraries(check_bitstring_distribution PRIVATE xacc::xacc tnqvm)