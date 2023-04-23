/*=================================================================
 * Implement the fast walsh-hadamard transform using C to accelerate the speed. 
 * fwht_user.c
 * Copyright 2020-2020 Haiwen Cao.
 *	
 *=================================================================*/
#include <stdio.h>
#include <string.h>
#include "mex.h"


/* The computational routine */
void fwht_user(double *u, mwSize N)
{
    mwSize i, j, k, ij;
    double temp;

    for(i=N>>1; i>0; i>>=1) {
        for(k=0; k<N; k += 2*i) {
            for(j=k, ij=i+k; j<k+i; j++, ij++) {
                temp   = u[j];
                u[j]  += u[ij];
                u[ij] = temp - u[ij];
            }
        }
    }

}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *inVector;               /* Nx1 input vector */
    size_t nrows;                   /* length of vector */
    double *outVector;              /* output vector */
    
    /* create a pointer to the real data in the input vector  */
    inVector = mxGetPr(prhs[0]);

    /* get length of the input vector */
    nrows = mxGetN(prhs[0]);
    
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)nrows, mxREAL);
    
  	/* get a pointer to the real data in the output vector  */
    outVector = mxGetPr(plhs[0]);
    
    memcpy(outVector, inVector, sizeof(double)*nrows);
    
    /* call the computational routine */
    fwht_user(outVector, (mwSize)nrows);
    

}
