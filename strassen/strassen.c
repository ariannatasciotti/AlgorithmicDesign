#include "matrix.h"

/*
 * this function performs the element-wise
 * subtraction of B from A and put the resulting
 * sub-matrix in C. The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub-matrix we want to deal with.
 */
void sub_matrix_blocks(float **C, float const *const *const A, float const *const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t rows, const size_t cols)
{
    for (size_t y = 0; y < rows; y++) {
      for (size_t x = 0; x < cols; x++) {
          C[y+C_f_row][x+C_f_col] =
               A[y+A_f_row][x+A_f_col] - B[y+B_f_row][x+B_f_col];
      }
    }
}

/*
 * this function performs the element-wise
 * sum of A and B and put the resulting
 * sub-matrix in C. The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub-matrix we want to deal with.
 */
void sum_matrix_blocks(float **C, float const *const *const A, float const *const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t rows, const size_t cols)
{
    for (size_t y = 0; y < rows; y++) {
      for (size_t x = 0; x < cols; x++) {
          C[y+C_f_row][x+C_f_col] =
               A[y+A_f_row][x+A_f_col] + B[y+B_f_row][x+B_f_col];
      }
    }
}


// row column product
void naive_aux(float **C, float const *const *const A, float const *const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t A_rows, const size_t A_cols, const size_t B_cols)
{
  for (size_t y = 0; y < A_rows; y++) {
    for (size_t x = 0; x < B_cols; x++) {
      float value = 0.0;
      for (size_t z = 0; z < A_cols; z++) {
        value += A[y + A_f_row][z + A_f_col]*B[z + B_f_row][x + B_f_col];
      }

      C[y + C_f_row][x + C_f_col] = value;
    }
  }
}



// Strassen algorithm with dynamic peeling
void strassen_aux(float **C, float const *const *const A, float const *const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t A_rows, const size_t A_cols, const size_t B_cols)
{
    if (A_rows < (1<<5) || A_cols < (1<<5) || B_cols < (1<<5))
    {
        naive_aux(C, A, B, C_f_row, C_f_col, A_f_row, A_f_col, B_f_row, B_f_col, A_rows, A_cols, B_cols);
        return;
    }
        
    size_t a_r0 = A_rows%2, a_c0 = A_cols%2, b_c0 = B_cols%2;
    
    size_t a_r1 = A_rows - a_r0, a_c1 = A_cols - a_c0, b_c1 = B_cols - b_c0;
      
    size_t a12_rows = a_r1, a12_cols = a_c0;
      
    size_t a21_rows = a_r0, a21_cols = a_c1;
      
    size_t a22_rows = a_r0, a22_cols = a_c0;
      
    size_t b12_rows = a_c1, b12_cols = b_c0;
      
    size_t b21_rows = a_c0, b21_cols = b_c1;
      
    size_t b22_rows = a_c0, b22_cols = b_c0;

    
    if(a_r0 || a_c0 || b_c0)
    {
        // call strassen on peeled matrices
        strassen_aux(C, A, B, C_f_row, C_f_col, A_f_row, A_f_col, B_f_row, B_f_col, a_r1, a_c1, b_c1);
        //printf("%d\t",0);
        
        // Fixedup
        
        // compute C21 = a21B11 + a22b21
        if (a_r0){  // a21B11 if A rows are odd
            
            naive_aux(C, A, B, C_f_row + a_r1, C_f_col, A_f_row + a_r1, A_f_col, B_f_row, B_f_col, a21_rows, a21_cols, b_c1);
        
            if (a_c0){  // also if A cols are odd
                
              float** C21 = allocate_matrix(a22_rows, b21_cols);
               
              // a22b21
              naive_aux(C21, A, B, 0, 0, A_f_row + a_r1, A_f_col + a_c1, B_f_row + a_c1, B_f_col, a22_rows, a22_cols, b21_cols);
                
              sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C21, C_f_row + a_r1, C_f_col, C_f_row + a_r1, C_f_col, 0, 0, a22_rows, b_c1);
                
              deallocate_matrix(C21, a22_rows);
            }
        }
        
        // compute C22 = a21b12 + a22b22
        if (a_r0 && b_c0){ // a21b12 if A rows and B cols are odd
            
            naive_aux(C, A, B, C_f_row + a_r1, C_f_col + b_c1, A_f_row + a_r1, A_f_col, B_f_row, B_f_col + b_c1, a21_rows, a21_cols, b12_cols);

            if (a_c0){ // also if A cols are odd
            
                float** C22 = allocate_matrix(a22_rows, b22_cols);
                
                // a22b22
                naive_aux(C22, A, B, 0, 0, A_f_row + a_r1, A_f_col + a_c1, B_f_row + a_c1, B_f_col + b_c1, a22_rows, a22_cols, b22_cols);
                
                sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C22, C_f_row + a_r1, C_f_col + b_c1, C_f_row + a_r1, C_f_col + b_c1, 0, 0, a22_rows, b22_cols);
                
                deallocate_matrix(C22, a22_rows);
            }
        }
            
        // compute C11 = A11B11 + a12b21
        if(a_c0){
            
          float** C11 = allocate_matrix(a12_rows, b21_cols);
            
          // a12b21
          naive_aux(C11, A, B, 0, 0, A_f_row, A_f_col + a_c1, B_f_row + a_c1, B_f_col, a12_rows, a12_cols, b21_cols);
            
          sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C11, C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, a_r1, b_c1);
            
          deallocate_matrix(C11, a12_rows);
        }
        
        // compute C12 = A11b12 + a12b22
        if (b_c0){ // A11b12 if B cols are odd
            
            naive_aux(C, A, B, C_f_row, C_f_col + b_c1, A_f_row, A_f_col, B_f_row, B_f_col + b_c1, a_r1, a_c1, b12_cols);
        
            if (a_c0){ // also if A cols are odd
                
                float** C12 = allocate_matrix(a12_rows, b22_cols);
                
                // a12b22
                naive_aux(C12, A, B, 0, 0, A_f_row, A_f_col + a_c1, B_f_row + a_c1, B_f_col + b_c1, a12_rows, a12_cols, b22_cols);
                
                sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C12, C_f_row, C_f_col + b_c1, C_f_row, C_f_col + b_c1, 0, 0, a12_rows, b22_cols);
                
                deallocate_matrix(C12, a12_rows);
            }
        }
        
            /*float** D = allocate_matrix(a12_rows, b21_cols);
            
            //a12*b21
            naive_aux(D, A, B, 0, 0, A_f_row, A_f_col + A_cols-1, B_f_row + A_cols-1, B_f_col, a12_rows, a12_cols, b21_cols);
            
            //C11 + a12*b21
            sum_matrix_blocks(C, D, C, C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, a12_rows, b21_cols);
            
            deallocate_matrix(D, a12_rows);
            
            // c12 =
            naive_aux(C, A, B, C_f_row, C_f_col + B_cols-1, 0, 0, 0, B_f_col + b_c1, a12_rows, A_cols, b12_cols);
        
            
            naive_aux(C, A, B,  C_f_row + A_rows-1, C_f_col, A_f_row + A_rows-1, A_f_col, B_f_row, B_f_col, a21_rows, A_cols, B_cols);*/
        
        return;
    }
    
    
    size_t A_r2 = A_rows/2;
    //printf("%d\t",A_r2);
    size_t A_c2 =  A_cols/2;
    //printf("%d\t",A_c2);
    size_t B_c2 = B_cols/2;
    //printf("%d\t",B_c2);
    
    

    float ***S_a = (float ***)malloc(sizeof(float **) * 5);
    for (size_t i = 0; i < 5; i++) {
        S_a[i] = allocate_matrix(A_r2, A_c2);
    }

    float ***S_b = (float ***)malloc(sizeof(float **) * 5);
    for (size_t i = 0; i < 5; i++) {
        S_b[i] = allocate_matrix(A_c2, B_c2);
    }

    float ***P = (float ***)malloc(sizeof(float **) * 7);
    for (size_t i = 0; i < 7; i++) {
        P[i] = allocate_matrix(A_r2, B_c2);
    }
    

    // S1 = B12 - B22
    sub_matrix_blocks(S_b[0], B, B, 0, 0, B_f_row, B_f_col + B_c2, B_f_row + A_c2, B_f_col + B_c2, A_c2, B_c2);

    // P1 = A11 x S1
    strassen_aux(P[0], A, (const float* const *const) S_b[0], 0, 0, A_f_row, A_f_col, 0, 0, A_r2,A_c2,B_c2);

    // S2 = A11 + A12
    sum_matrix_blocks(S_a[0], A, A, 0, 0, A_f_row, A_f_col, A_f_row, A_f_col + A_c2, A_r2,A_c2);
    
    // P2 = S2 x B22
    strassen_aux(P[1], (const float* const *const) S_a[0], B, 0, 0, 0, 0, B_f_row + A_c2, B_f_col + B_c2, A_r2, A_c2, B_c2);

    // S3 = A21 + A22
    sum_matrix_blocks(S_a[1], (const float* const *const) A, A, 0, 0, A_f_row + A_r2, A_f_col, A_f_row + A_r2, A_f_col + A_c2, A_r2, A_c2);
    
    // P3 = S3 x B11
    strassen_aux(P[2], (const float* const *const) S_a[1], B, 0, 0, 0, 0, B_f_row, B_f_col, A_r2, A_c2, B_c2);
    
    // S4 = B21 - B11
    sub_matrix_blocks(S_b[1], (const float* const *const) B, B, 0, 0, B_f_row + A_c2, B_f_col, B_f_row, B_f_col, A_c2, B_c2);
    
    // P4 = A22 x S4
    strassen_aux(P[3], A, (const float* const *const) S_b[1], 0, 0, A_f_row + A_r2, A_f_col + A_c2, 0, 0, A_r2, A_c2, B_c2);

    // S5 = A11 + A22
    sum_matrix_blocks(S_a[2], (const float* const *const) A, A, 0, 0, A_f_row, A_f_col, A_f_row + A_r2, A_f_col + A_c2, A_r2, A_c2);

    // S6 = B11 + B22
    sum_matrix_blocks(S_b[2], (const float* const *const) B, B, 0, 0, B_f_row, B_f_col, B_f_row + A_c2, B_f_col + B_c2, A_c2, B_c2);

    // P5 = S5 x S6
    strassen_aux(P[4], (const float* const *const) S_a[2], (const float* const *const) S_b[2], 0, 0, 0, 0, 0, 0, A_r2, A_c2, B_c2);

    // S7 = A12 - A22
    sub_matrix_blocks(S_a[3], (const float* const *const) A, A, 0, 0, A_f_row, A_f_col + A_c2, A_f_row + A_r2, A_f_col + A_c2, A_r2, A_c2);
    
    // S8 = B21 + B22
    sum_matrix_blocks(S_b[3], (const float* const *const) B, B, 0, 0, B_f_row + A_c2, B_f_col, B_f_row + A_c2, B_f_col + B_c2, A_c2, B_c2);

    // P6 = S7 x S8
    strassen_aux(P[5], (const float* const *const) S_a[3], (const float* const *const) S_b[3], 0, 0, 0, 0, 0, 0, A_r2, A_c2, B_c2);

    // S9 = A11 - A21
    sub_matrix_blocks(S_a[4], (const float* const *const) A, A, 0, 0, A_f_row, A_f_col, A_f_row + A_r2, A_f_col, A_r2, A_c2);
    
    // S10 = B11 + B12
    sum_matrix_blocks(S_b[4], (const float* const *const) B, B, 0, 0, B_f_row, B_f_col, B_f_row, B_f_col + B_c2, A_c2, B_c2);

    // P7 = S9 x S10
    strassen_aux(P[6], (const float* const *const) S_a[4], (const float* const *const) S_b[4], 0, 0, 0, 0, 0, 0, A_r2, A_c2, B_c2);

    // C11 = P5 + P4 - P2 + P6
    sum_matrix_blocks(C, (const float* const *const) P[4], (const float* const *const) P[3], C_f_row, C_f_col,0, 0, 0, 0, A_r2, B_c2);
    
    sub_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[1], C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, A_r2, B_c2);
    
    sum_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[5], C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, A_r2, B_c2);

    // C12 = P1 + P2
    sum_matrix_blocks(C, (const float* const *const) P[0], (const float* const *const) P[1], C_f_row, C_f_col+B_c2, 0, 0, 0, 0, A_r2, B_c2);

    // C21 = P3 + P4
    sum_matrix_blocks(C, (const float* const *const) P[2], (const float* const *const) P[3], C_f_row+A_r2, C_f_col, 0, 0, 0, 0, A_r2, B_c2);

    // C22 = P5 + P1 - P3 - P7
    sum_matrix_blocks(C, (const float* const *const) P[4], (const float* const *const) P[0], C_f_row+A_r2, C_f_col+B_c2, 0, 0, 0, 0, A_r2, B_c2);
    
    sub_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[2], C_f_row+A_r2, C_f_col+B_c2, C_f_row+A_r2, C_f_col+B_c2, 0, 0, A_r2, B_c2);
    
    sub_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[6], C_f_row+A_r2, C_f_col+B_c2, C_f_row+A_r2, C_f_col+B_c2, 0, 0, A_r2, B_c2);
    
    
    for (size_t i = 0; i < 5; i++) {
      deallocate_matrix(S_a[i], A_r2);
    }
    free(S_a);
    
    for (size_t i = 0; i < 5; i++) {
      deallocate_matrix(S_b[i], A_c2);
    }
    free(S_b);
    
    for (size_t i = 0; i < 7; i++) {
      deallocate_matrix(P[i], A_r2);
    }
    free(P);
}


void strassen_matrix_multiplication(float **C, float const *const *const A, float const *const *const B, size_t A_rows, size_t A_cols, size_t B_cols)
{

  strassen_aux(C, A, B, 0, 0, 0, 0, 0, 0, A_rows, A_cols, B_cols);

}


// Strassen algorithm with dynamic peeling and reduction of memory allocations
void strassen_aux_imp(float **C, float const *const *const A, float const *const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t A_rows, const size_t A_cols, const size_t B_cols)
{
    if (A_rows < (1<<5) || A_cols < (1<<5) || B_cols < (1<<5))
    {
        naive_aux(C, A, B, C_f_row, C_f_col, A_f_row, A_f_col, B_f_row, B_f_col, A_rows, A_cols, B_cols);
        return;
    }
    
    size_t a_r0 = A_rows%2, a_c0 = A_cols%2, b_c0 = B_cols%2;
    
    size_t a_r1 = A_rows - a_r0, a_c1 = A_cols - a_c0, b_c1 = B_cols - b_c0;
      
    size_t a12_rows = a_r1, a12_cols = a_c0;
      
    size_t a21_rows = a_r0, a21_cols = a_c1;
      
    size_t a22_rows = a_r0, a22_cols = a_c0;
      
    size_t b12_rows = a_c1, b12_cols = b_c0;
      
    size_t b21_rows = a_c0, b21_cols = b_c1;
      
    size_t b22_rows = a_c0, b22_cols = b_c0;

    
    if(a_r0 || a_c0 || b_c0)
    {
        // call strassen on peeled matrices
        strassen_aux_imp(C, A, B, C_f_row, C_f_col, A_f_row, A_f_col, B_f_row, B_f_col, a_r1, a_c1, b_c1);
        //printf("%d\t",0);
        
        // Fixedup
        
        // compute C21 = a21B11 + a22b21
        if (a_r0){  // a21B11 if A rows are odd
            
            naive_aux(C, A, B, C_f_row + a_r1, C_f_col, A_f_row + a_r1, A_f_col, B_f_row, B_f_col, a21_rows, a21_cols, b_c1);
        
            if (a_c0){  // also if A cols are odd
                
              float** C21 = allocate_matrix(a22_rows, b21_cols);
               
              // a22b21
              naive_aux(C21, A, B, 0, 0, A_f_row + a_r1, A_f_col + a_c1, B_f_row + a_c1, B_f_col, a22_rows, a22_cols, b21_cols);
                
              sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C21, C_f_row + a_r1, C_f_col, C_f_row + a_r1, C_f_col, 0, 0, a22_rows, b_c1);
                
              deallocate_matrix(C21, a22_rows);
            }
        }
        
        // compute C22 = a21b12 + a22b22
        if (a_r0 && b_c0){ // a21b12 if A rows and B cols are odd
            
            naive_aux(C, A, B, C_f_row + a_r1, C_f_col + b_c1, A_f_row + a_r1, A_f_col, B_f_row, B_f_col + b_c1, a21_rows, a21_cols, b12_cols);

            if (a_c0){ // also if A cols are odd
            
                float** C22 = allocate_matrix(a22_rows, b22_cols);
                
                // a22b22
                naive_aux(C22, A, B, 0, 0, A_f_row + a_r1, A_f_col + a_c1, B_f_row + a_c1, B_f_col + b_c1, a22_rows, a22_cols, b22_cols);
                
                sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C22, C_f_row + a_r1, C_f_col + b_c1, C_f_row + a_r1, C_f_col + b_c1, 0, 0, a22_rows, b22_cols);
                
                deallocate_matrix(C22, a22_rows);
            }
        }
            
        // compute C11 = A11B11 + a12b21
        if(a_c0){
            
          float** C11 = allocate_matrix(a12_rows, b21_cols);
            
          // a12b21
          naive_aux(C11, A, B, 0, 0, A_f_row, A_f_col + a_c1, B_f_row + a_c1, B_f_col, a12_rows, a12_cols, b21_cols);
            
          sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C11, C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, a_r1, b_c1);
            
          deallocate_matrix(C11, a12_rows);
        }
        
        // compute C12 = A11b12 + a12b22
        if (b_c0){ // A11b12 if B cols are odd
            
            naive_aux(C, A, B, C_f_row, C_f_col + b_c1, A_f_row, A_f_col, B_f_row, B_f_col + b_c1, a_r1, a_c1, b12_cols);
        
            if (a_c0){ // also if A cols are odd
                
                float** C12 = allocate_matrix(a12_rows, b22_cols);
                
                // a12b22
                naive_aux(C12, A, B, 0, 0, A_f_row, A_f_col + a_c1, B_f_row + a_c1, B_f_col + b_c1, a12_rows, a12_cols, b22_cols);
                
                sum_matrix_blocks(C, (const float* const* const)C, (const float* const* const)C12, C_f_row, C_f_col + b_c1, C_f_row, C_f_col + b_c1, 0, 0, a12_rows, b22_cols);
                
                deallocate_matrix(C12, a12_rows);
            }
        }
        
        return;
    }
    
    
    size_t A_r2 = A_rows/2;
    //printf("%d\t",A_r2);
    size_t A_c2 =  A_cols/2;
    //printf("%d\t",A_c2);
    size_t B_c2 = B_cols/2;
    //printf("%d\t",B_c2);
    

    float** S_a = allocate_matrix(A_r2, A_c2);

    float** S_b = allocate_matrix(A_c2, B_c2);

    float ***P = (float ***)malloc(sizeof(float **) * 3);
    for (size_t i = 0; i < 3; i++) {
        P[i] = allocate_matrix(A_r2, B_c2);
    }


    // S1 = B12 - B22
    sub_matrix_blocks(S_b, B, B, 0, 0, B_f_row, B_f_col + B_c2, B_f_row + A_c2, B_f_col + B_c2, A_c2, B_c2);

    // P1 = A11 x S1
    strassen_aux_imp(P[0], A, (const float* const *const) S_b, 0, 0, A_f_row, A_f_col, 0, 0, A_r2,A_c2,B_c2);

    // S2 = A11 + A12
    sum_matrix_blocks(S_a, A, A, 0, 0, A_f_row, A_f_col, A_f_row, A_f_col + A_c2, A_r2,A_c2);

    // P2 = S2 x B22
    strassen_aux_imp(P[1], (const float* const *const) S_a, B, 0, 0, 0, 0, B_f_row + A_c2, B_f_col + B_c2, A_r2, A_c2, B_c2);

    // C12 = P1 + P2
    sum_matrix_blocks(C, (const float* const *const) P[0], (const float* const *const) P[1], C_f_row, C_f_col+B_c2, 0, 0, 0, 0, A_r2, B_c2);
        
    // S3 = A21 + A22
    sum_matrix_blocks(S_a, A, A, 0, 0, A_f_row + A_r2, A_f_col, A_f_row + A_r2, A_f_col + A_c2, A_r2, A_c2);

    // P3 = S3 x B11
    strassen_aux_imp(P[2], (const float* const *const) S_a, B, 0, 0, 0, 0, B_f_row, B_f_col, A_r2, A_c2, B_c2);

    sub_matrix_blocks(C, (const float* const *const) P[0], (const float* const *const) P[2], C_f_row+A_r2, C_f_col+B_c2, 0, 0, 0, 0, A_r2, B_c2);
        
    // S4 = B21 - B11
    sub_matrix_blocks(S_b, B, B, 0, 0, B_f_row + A_c2, B_f_col, B_f_row, B_f_col, A_c2, B_c2);

    // P4 = A22 x S4
    strassen_aux_imp(P[0], A, (const float* const *const) S_b, 0, 0, A_f_row + A_r2, A_f_col + A_c2, 0, 0, A_r2, A_c2, B_c2);
        
    // C21 = P3 + P4
    sum_matrix_blocks(C, (const float* const *const) P[0], (const float* const *const) P[2], C_f_row+A_r2, C_f_col, 0, 0, 0, 0, A_r2, B_c2);

    // C11 = P5 + P4 - P2 + P6
    sub_matrix_blocks(C, (const float* const *const) P[0], (const float* const *const) P[1], C_f_row, C_f_col, 0, 0, 0, 0, A_r2, B_c2);
    
    // S5 = A11 + A22
    sum_matrix_blocks(S_a, A, A, 0, 0, A_f_row, A_f_col, A_f_row + A_r2, A_f_col + A_c2, A_r2, A_c2);

    // S6 = B11 + B22
    sum_matrix_blocks(S_b, B, B, 0, 0, B_f_row, B_f_col, B_f_row + A_c2, B_f_col + B_c2, A_c2, B_c2);

    // P5 = S5 x S6
    strassen_aux_imp(P[0], (const float* const *const) S_a, (const float* const *const) S_b, 0, 0, 0, 0, 0, 0, A_r2, A_c2, B_c2);
        
    sum_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[0], C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, A_r2, B_c2);
    
    sum_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[0], C_f_row+A_r2, C_f_col+B_c2, C_f_row+A_r2, C_f_col+B_c2, 0, 0, A_r2, B_c2);

    // S7 = A12 - A22
    sub_matrix_blocks(S_a, A, A, 0, 0, A_f_row, A_f_col + A_c2, A_f_row + A_r2, A_f_col + A_c2, A_r2, A_c2);

    // S8 = B21 + B22
    sum_matrix_blocks(S_b, B, B, 0, 0, B_f_row + A_c2, B_f_col, B_f_row + A_c2, B_f_col + B_c2, A_c2, B_c2);

    // P6 = S7 x S8
    strassen_aux_imp(P[0], (const float* const *const) S_a, (const float* const *const) S_b, 0, 0, 0, 0, 0, 0, A_r2, A_c2, B_c2);
    
    sum_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[0], C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, A_r2, B_c2);

    // S9 = A11 - A21
    sub_matrix_blocks(S_a, A, A, 0, 0, A_f_row, A_f_col, A_f_row + A_r2, A_f_col, A_r2, A_c2);

    // S10 = B11 + B12
    sum_matrix_blocks(S_b, B, B, 0, 0, B_f_row, B_f_col, B_f_row, B_f_col + B_c2, A_c2, B_c2);

    // P7 = S9 x S10
    strassen_aux_imp(P[0], (const float* const *const) S_a, (const float* const *const) S_b, 0, 0, 0, 0, 0, 0, A_r2, A_c2, B_c2);
    
    sub_matrix_blocks(C, (const float* const *const) C, (const float* const *const) P[0], C_f_row+A_r2, C_f_col+B_c2, C_f_row+A_r2, C_f_col+B_c2, 0, 0, A_r2, B_c2);
        
        
    deallocate_matrix(S_a, A_r2);

    deallocate_matrix(S_b, A_c2);

    for (size_t i = 0; i < 3; i++) {
      deallocate_matrix(P[i], A_r2);
    }
    free(P);
}
    
void strassen_matrix_multiplication_imp(float **C, float const *const *const A, float const *const *const B, size_t A_rows, size_t A_cols, size_t B_cols)
{

  strassen_aux_imp(C, A, B, 0, 0, 0, 0, 0, 0, A_rows, A_cols, B_cols);

}
