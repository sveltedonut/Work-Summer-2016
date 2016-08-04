//
//  inprod_part0_6.c
//  
//
//  Created by Harry Jiang on 2016-05-12.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define S_OUT       plhs[0]
    
    // Macros for input
    #define S_IN        prhs[0]
    #define WIDTH       prhs[1]
    #define FX1         prhs[2]
    #define FX2         prhs[3]
    #define TNM         prhs[4]
    #define ITER        prhs[5]
    
    double *s_out, *tnm, *fx1, *fx2, *width, *iter;
    int *D, i, j, K;
    
    // get needed inputs
    tnm     = mxGetPr(TNM);
    width   = mxGetPr(WIDTH);
    fx1     = mxGetPr(FX1);
    fx2     = mxGetPr(FX2);
    iter    = mxGetPr(ITER);
    
    // set up s
    S_OUT   = S_IN;
    s_out   = mxGetPr(S_OUT);
    D       = mxGetDimensions(S_IN);
    K       = mxGetM(FX1);
    
    // assign (iter)th set of solutions
    for (i = 0; i < D[1]; i++) {
        for (j = 0; j < D[2]; j++) {
            int k, index;
            double n = 0;
            // multiply fx1^T and fx2
            for (k = 0; k < K; k++){
                n += fx1[k + K*i] * fx2[k + K*j];
            }
            index = (int)iter[0] - 1 + D[0]*i + D[0]*D[1]*j;
            s_out[index] = (s_out[index - 1] + (n * width[0]/tnm[0])) / 2;
        }
    }
    
    return;
}