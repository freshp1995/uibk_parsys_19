//from https://rosettacode.org/wiki/N-queens_problem#C

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define MAXN 31

void calc(int n, int limit, int *num, int *cols, int *diagl, int *diagr, int *posibs) {
	for (int q0=0; q0<limit -2; q0++) {
			for (int q1=q0+2; q1<n; q1++){
			int bit0 = 1<<q0;
			int bit1 = 1<<q1;

			int d=0; // d is our depth in the backtrack stack 
			cols[0] = bit0 | bit1 | (-1<<n); // The -1 here is used to fill all 'coloumn' bits after n ...
			diagl[0]= (bit0<<1 | bit1)<<1;
			diagr[0]= (bit0>>1 | bit1)>>1;
		
			//  The variable posib contains the bitmask of possibilities we still have to try in a given row ...
			int posib = ~(cols[0] | diagl[0] | diagr[0]);
		
			while (d >= 0) {
				while(posib) {
				int bit = posib & -posib; // The standard trick for getting the rightmost bit in the mask
				int ncols= cols[d] | bit;
				int ndiagl = (diagl[d] | bit) << 1;
				int ndiagr = (diagr[d] | bit) >> 1;
				int nposib = ~(ncols | ndiagl | ndiagr);
				posib^=bit; // Eliminate the tried possibility.
		
				// The following is the main additional trick here, as recognizing solution can not be done using stack level (d),
				// since we save the depth+backtrack time at the end of the enumeration loop. However by noticing all coloumns are
				// filled (comparison to -1) we know a solution was reached ...
				// Notice also that avoiding an if on the ncols==-1 comparison is more efficient!
				(*num) += ncols==-1; 
				
		
				if (nposib) {
					if (posib) { // This if saves stack depth + backtrack operations when we passed the last possibility in a row.
					posibs[d++] = posib; // Go lower in stack ..
					}
					cols[d] = ncols;
					diagl[d] = ndiagl;
					diagr[d] = ndiagr;
					posib = nposib;
				}
				}
				posib = posibs[--d]; // backtrack ...
			}
			}
		}

}
 
int nqueens(int n)
{
	int cols[MAXN], diagl[MAXN], diagr[MAXN], posibs[MAXN]; // Our backtracking 'stack' 
	int num=0;
	//
	// The top level is two fors, to save one bit of symmetry in the enumeration by forcing second queen to
	// be AFTER the first queen.
	//

	int num_threads = omp_get_num_threads();
	int balance = n / num_threads;

	
		

	//#pragma omp parallel for shared(num, cols, diagl, diagr, posibs) firstprivate(n, balance), default(none)
	for (int i = 0; i < n; i+=balance) {

		int limit = i + balance > n ? n : i + balance;
		//#pragma omp task shared(num, cols, diagl, diagr, posibs) firstprivate(n, limit) default(none)
		calc(n, limit, &num, cols, diagl, diagr, posibs);
	}
		

	
	return num*2;
}
 
 
int main(int argc, char **argv)
{
  if(argc != 2) {
    printf("usage: number of queens n\n");
    return 1;
  }
  
  int n = atoi(argv[1]);
  if(n<1 || n > MAXN) {
    printf("n must be between 2 and 31!\n");
  }
  
  clock_t begin = clock();
  int result = nqueens(n);
  clock_t end = clock();

  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("%f\n", time_spent);
  //printf("Number of solution for %d is %d\n",n,result);
  
  return EXIT_SUCCESS;
  
}
 
