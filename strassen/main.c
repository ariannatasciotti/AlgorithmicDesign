#include <stdio.h>

#include "test.h"
#include "matrix.h"
#include "strassen.h"


int main(int argc, char *argv[])
{
  
  size_t n = (1 << 12)+3;

  float **A = allocate_random_matrix(n, n);
    
  float **B = allocate_random_matrix(n, n);
    
  float **C0 = allocate_matrix(n, n);
    
  float **C1 = allocate_matrix(n, n);
    
  float **C2 = allocate_matrix(n, n);
    
  //printf("%d\t\tn_rows",n);
  
  printf("Strassen's Alg.\tStrassen Imp\tNaive Alg.\tSame results\tsize\n");
    
  for (size_t j = 2; j <= n; j *= 2)
  {
    
    // generalized strassen for non square and odd matrices
    printf("%lf\t", test(strassen_matrix_multiplication, C1, A, B, j+3, j, j-1));
    fflush(stdout);
      
    // generalized strassen for non square and odd matrices and reduction of memory allocations
    printf("%lf\t", test(strassen_matrix_multiplication_imp, C2, A, B, j+3, j,j-1));
    fflush(stdout);
      
    // row column multiplication
    printf("%lf\t", test(naive_matrix_multiplication, C0, A, B, j+3, j, j-1));
    fflush(stdout);
      
    // check if the results are correct: print 1 only if strassen alg and strassen imp give the same results obtained with naive
    printf("%d\t\t ", same_matrix((float const *const *const)C1, (float const *const *const)C2, j+3, j-1) && same_matrix((float const *const *const)C0, (float const *const *const)C2, j+3, j-1) && same_matrix((float const *const *const)C1, (float const *const *const)C0, j+3, j-1));
    fflush(stdout);
      
    // size of result matrix
    printf("%ld\n", (j+3)*(j-1));
    fflush(stdout);
  }

  deallocate_matrix(A, n);
    
  deallocate_matrix(B, n);
    
  deallocate_matrix(C0, n);
    
  deallocate_matrix(C1, n);
    
  deallocate_matrix(C2, n);

  return 0;
}
