project(svd)

add_library(svd "")

set_target_properties(svd PROPERTIES LINKER_LANGUAGE CXX)

target_sources(
    svd
    PRIVATE
        Singular_Value_Decomposition_Preamble.hpp 
        Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp 
        Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp 
        Singular_Value_Decomposition_Main_Kernel_Body.hpp
    PUBLIC
)