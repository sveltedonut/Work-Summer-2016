//
//  inprod_part0_5.c
//  
//
//  Created by Harry Jiang on 2016-05-12.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define TNM_OUT         plhs[0]
    #define X_OUT           plhs[1]
    
    // Macros for input
    #define TNM_IN          prhs[0]
    #define RNGI            prhs[1]
    #define WIDTH           prhs[2]
    
    double *width, *tnm, delta, *x, *rngi;
    int i;
    
    TNM_OUT = TNM_IN;
    tnm     = mxGetPr(TNM_OUT);
    tnm[0] *= 2;
    
    width   = mxGetPr(WIDTH);
    delta   = width[0]/tnm[0];
    
    X_OUT   = mxCreateDoubleMatrix(1, (int)tnm[0], mxREAL);
    x       = mxGetPr(X_OUT);
    
    rngi     = mxGetPr(RNGI);
    
    x[0]    = rngi[0] + delta/2;
    for (i = 1; i < (int)tnm[0]; i++){
        x[i] = x[i-1] + delta;
    }
    
    return;
}
