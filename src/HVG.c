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
void sumSeq(int *start, int *size, int *sumVect)
{
    /*
    This function provides a simple sequential sum
    where F[n] = F[n-1] + n
    */
    int i, j ;
    j = 0 ;
    for(i = *start; i < (*start + *size); i++){
        if(i == *start){
            sumVect[j] = i ;
        }
        else{
            sumVect[j] = sumVect[j-1] + i ;
        }
        j ++ ;
    }
}

/* This is just cheat-sheet from C: http://mcglinn.web.unc.edu/blog/linking-c-with-r-in-windows/ */
void fiboSeq(int *size, int *sumVect)
{
    /*
    This function returns the Fibonacci sequence
    where F[n] = F[n-1] + F[n-2]
    */
    int i ;
    sumVect[0] = 0 ;
    sumVect[1] = 1 ;
    for(i = 2; i < *size; i++){
        sumVect[i] = sumVect[i-1] + sumVect[i-2] ;
    }
}

/* very simple  */
void HVG_C(double *x, int *N, int *forwardBlocker, int *backwardBlocker)
{

    int i, j;

    /* look forward to first blocker, then stop: */
    for (i = 0; i < (*N-1); i++) {
        for(j = i + 1; j < (*N); j++) { // equivalent to which statement in R
            if(x[j] >= x[i]) {
                forwardBlocker[i] = j + 1; // +1 is to correct for R indexing
                break;
            }
        }
    }

    /* look backward to the first hit, then stop: */
    for (i = (*N-1); i > 0; i--) {
        for(j = i - 1; j > -1; j--) { // equivalent to which statement in R
            if(x[i] <= x[j]) {
                backwardBlocker[i] = j +1; // +1 is to correct for R indexing
                break;
            }
        }
    }
}


