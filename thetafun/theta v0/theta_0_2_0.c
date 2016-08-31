//
//  theta_0_2_0.c
//  
//
//  Created by Harry Jiang on 2016-08-15.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define STEPITER_OUT    plhs[0]
    #define TINANEW_OUT     plhs[1]
    
    // Macros for input
    #define STEPITER_IN prhs[0]
    #define STEP        prhs[1]
    #define SCRMIN      prhs[2]
    #define SCRMAX      prhs[3]
    #define FAILINDEX   prhs[4]
    #define FAC         prhs[5]
    #define TINANEW_IN  prhs[6]
    #define ACTIVE      prhs[7]
    
    long M          = mxGetM(STEP);
    
    double *step    = mxGetPr(STEP);
    double stepnew;
    
    STEPITER_OUT    = STEPITER_IN;
    double *siter   = mxGetPr(STEPITER_OUT);
    siter[0]       += 1;
    
    double *active  = mxGetPr(ACTIVE);
    
    TINANEW_OUT     = TINANEW_IN;
    double *tinanew = mxGetPr(TINANEW_OUT);
    double tinanewfail;
    
    double scrmin   = mxGetPr(SCRMIN)[0];
    double scrmax   = mxGetPr(SCRMAX)[0];
    
    double fac      = mxGetPr(FAC)[0];
    
    mxLogical *fail = mxGetLogicals(FAILINDEX);
    
    long f = 0;
    long m;
    for (m = 0; m < M; m++) {
        // Compute current step size
        if (fail[m]) {
            stepnew     = fac * step[m];
            tinanewfail = tinanew[m];
            if (tinanewfail - stepnew <= scrmin) stepnew = tinanewfail - scrmin * 1.01;
            if (tinanewfail - stepnew >= scrmax) stepnew = -(scrmax * 0.99 - tinanewfail);
            
            // Update ability values
            tinanew[m] -= stepnew;
            
            f++;
        }
    }
    
    return;
}