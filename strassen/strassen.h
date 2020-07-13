#ifndef __STRASSEN__

void strassen_matrix_multiplication(float **C, float const *const *const A, float const *const *const B, const size_t A_rows, const size_t A_cols, const size_t B_cols);

void strassen_matrix_multiplication_imp(float **C, float const *const *const A, float const *const *const B, size_t A_rows, size_t A_cols, size_t B_cols);

#endif //__STRASSEN__
