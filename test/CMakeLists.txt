add_executable(
    unit_tests
    rasterization_tests.cpp
)

target_link_libraries(
    unit_tests
    gtest_main
    solverlib
)

add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/unit_tests
)