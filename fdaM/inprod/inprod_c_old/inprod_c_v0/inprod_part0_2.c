//
//  inprod_part0_2.c
//  
//
//  Created by Harry Jiang on 2016-05-10.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define H_OUT       plhs[0]
    #define S_OUT       plhs[1]
    
    // Macros for input
    #define JMAXP       prhs[0]
    #define NREP1       prhs[1]
    #define NREP2       prhs[2]
    
    double *h, *s, *nrep1, *nrep2;
    int *jmaxp, j;
    
    // retrieve JMAXP, nrep1, nrep2
    jmaxp   = (int *)mxGetPr(JMAXP);
    nrep1   = mxGetPr(NREP1);
    nrep2   = mxGetPr(NREP2);
    
    // create uninitialized h array
    H_OUT = mxCreateDoubleMatrix(jmaxp[0], 1, mxREAL);
    //mxSetM(H_OUT, jmaxp_i);
    //mxSetN(H_OUT, 1);
    //mxSetData(H_OUT, mxMalloc(sizeof(double)*jmaxp_i));
    h = mxGetPr(H_OUT);
    
    // assign values to h
    h[0] = 1;
    h[1] = 0.25;
    for (j = 2; j < jmaxp[0]; j++) {
        h[j] = 1;
    }
    
    mwSize dims[3]  = {jmaxp[0], (int)nrep1[0], (int)nrep2[0]};
    S_OUT           = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    
    return;
}