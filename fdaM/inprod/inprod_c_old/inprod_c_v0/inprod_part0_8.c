//
//  inprod_part0_8.c
//  
//
//  Created by Harry Jiang on 2016-05-13.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define S_OUT       plhs[0]
    #define H_OUT       plhs[1]
    
    // Macros for input
    #define S_IN        prhs[0]
    #define H_IN        prhs[1]
    #define ITER        prhs[2]
    
    double *s, *h;
    int iter, *D, j, k, index;
    
    iter = (int)mxGetPr(ITER)[0];
    
    S_OUT = S_IN;
    H_OUT = H_IN;
    s     = mxGetPr(S_OUT);
    h     = mxGetPr(H_OUT);
    D    = mxGetDimensions(S_IN);
    
    for (j = 0; j < D[1]; j++) {
        for (k = 0; k < D[2]; k++) {
            index = D[0]*j + D[0]*D[1]*k;
            s[iter + index] = s[iter + index - 1];
        }
    }
    
    h[iter] = 0.25 * h[iter - 1];
    
    return;
}