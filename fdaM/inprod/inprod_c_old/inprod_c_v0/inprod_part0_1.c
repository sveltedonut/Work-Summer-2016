//
//  inprod_part0_1.c
//  
//
//  Created by Harry Jiang on 2016-05-10.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define RNGI        plhs[0]
    #define ITER        plhs[1]
    #define WIDTH       plhs[2]
    #define JMAXP       plhs[3]
    
    // Macros for input
    #define NRNG        prhs[0]
    #define IRNG        prhs[1]
    #define RNGVEC      prhs[2]
    #define JMAX        prhs[3]
    
    double *rngi, *width, *nrng, *irng, *rngvec, *jmax;
    int *iter, *jmaxp, i, n;
    
    // retrieve irng and nrng
    irng    = mxGetPr(IRNG);
    nrng    = mxGetPr(NRNG);
    i       = (int)irng[0];
    n       = (int)nrng[0];
    
    // retrieve rngvec
    rngvec  = mxGetPr(RNGVEC);
    
    // initialize and assign rngi
    RNGI    = mxCreateDoubleMatrix(1, 2, mxREAL);
    rngi    = mxGetPr(RNGI);
    rngi[0] = rngvec[i-2];
    rngi[1] = rngvec[i-1];
    
    // adjust rngi
    if (i > 2) {
        rngi[0] += 1e-10;
    }
    if (i < n){
        rngi[1] -= 1e-10;
    }
    
    // initialize and assign iter
    ITER    = mxCreateDoubleMatrix(1, 1, mxREAL);
    iter    = (int *)mxGetPr(ITER);
    iter[0] = 1;
    
    // initialize ans calculate width
    WIDTH       = mxCreateDoubleMatrix(1, 1, mxREAL);
    width       = mxGetPr(WIDTH);
    width[0]    = rngi[1] - rngi[0];
    
    // initialize and calculate JMAXP
    jmax        = mxGetPr(JMAX);
    JMAXP       = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    jmaxp       = (int *)mxGetPr(JMAXP);
    jmaxp[0]    = (int)jmax[0] + 1;
    
    return;
}