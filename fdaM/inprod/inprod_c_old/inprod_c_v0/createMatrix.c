//
//  createMatrix.c
//  
//
//  Created by Harry Jiang on 2016-05-10.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define OUT             plhs[0]
    
    // Macros for input
    #define INA             prhs[0]
    #define INB             prhs[1]

    double *A, *B, *C;
    int M, i;
    
    A = mxGetPr(INA);
    B = mxGetPr(INB);
    M = mxGetM(INA);
    OUT = mxCreateDoubleMatrix(M, 1, mxREAL);
    C = mxGetPr(OUT);
    
    for (i = 0; i < M; i++){
        C[i] = A[i] + B[i];
    }
    
    return;
}
