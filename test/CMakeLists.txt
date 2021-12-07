# Parent project does not export its library target, so this CML implicitly
# depends on being added from it, i.e. the testing is done only from the build
# tree and is not feasible from an install location

project(atomic-potentialsTests LANGUAGES CXX)

add_executable(atomic-potentials_test source/atomic-potentials_test.cpp)
target_link_libraries(atomic-potentials_test PRIVATE atomic-potentials_lib)
target_compile_features(atomic-potentials_test PRIVATE cxx_std_17)

add_test(NAME atomic-potentials_test COMMAND atomic-potentials_test)
