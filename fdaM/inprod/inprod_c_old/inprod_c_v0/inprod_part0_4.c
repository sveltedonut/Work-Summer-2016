//
//  inprod_part0_4.c
//  
//
//  Created by Harry Jiang on 2016-05-11.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define S_OUT       plhs[0]
    #define TNM         plhs[1]
    
    // Macros for input
    #define S_IN        prhs[0]
    #define WIDTH       prhs[1]
    #define FX1         prhs[2]
    #define FX2         prhs[3]
    
    double *s_out, *tnm, *fx1, *fx2, *width;
    int *D, i, j, K;
    
    // initialize tnm
    TNM     = mxCreateDoubleMatrix(1, 1, mxREAL);
    tnm     = mxGetPr(TNM);
    tnm[0]  = 0.5;
    
    // get needed inputs
    width   = mxGetPr(WIDTH);
    fx1     = mxGetPr(FX1);
    fx2     = mxGetPr(FX2);
    
    // set up s
    S_OUT   = S_IN;
    s_out   = mxGetPr(S_OUT);
    D       = mxGetDimensions(S_IN);
    K       = mxGetM(FX1);
    //mexPrintf("N(FX1) = %d\n", mxGetN(FX1));
    //mexPrintf("N(FX2) = %d\n", mxGetN(FX2));
    
    // assign first set of solutions
    for (i = 0; i < D[1]; i++) {
        for (j = 0; j < D[2]; j++) {
            int k;
            double n = 0;
            // multiply fx1^T and fx2
            for (k = 0; k < K; k++){
                n += fx1[k + K*i] * fx2[k + K*j];
            }
            s_out[D[0]*i + D[0]*D[1]*j] = n*width[0]/2;
        }
    }
    
    return;
}