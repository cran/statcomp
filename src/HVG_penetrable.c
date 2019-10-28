#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rdefines.h>
/*
Filename: "sequence_examples.c"
Return a vectors of sequentially summed values
Arguments:
start -- value to start the sum at
size -- the number of elements to return
sumVect -- the vector of summed output values
*/
void HVG_penetrable_C(double *x, int *N, int *forwardBlocker, int *backwardBlocker, int *rho){

    int i, j, cur_rho;
    cur_rho = 0;

    /* look forward to first blocker, then stop: */
    for (i = 0; i < (*N-1); i++) {
            cur_rho = 0;
        for(j = i + 1; j < (*N); j++) { // equivalent to which statement in R
            if(x[j] >= x[i]) {
                forwardBlocker[cur_rho * (*N) + i] = j + 1; // +1 is to correct for R indexing
                cur_rho = cur_rho + 1;
                // break;
                if(cur_rho > (*rho)) {
                    break;
                }
            }
        }
    }

    /* look backward to the first hit, then stop: */
    for (i = (*N-1); i > 0; i--) {
            cur_rho = 0;
        for(j = i - 1; j > -1; j--) { // equivalent to which statement in R
            if(x[i] <= x[j]) {
                backwardBlocker[cur_rho * (*N) + i] = j +1; // +1 is to correct for R indexing
                cur_rho = cur_rho + 1;
                // break;
                if(cur_rho > *rho) {
                    break;
                }
            }
        }
    }
}

