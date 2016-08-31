//
//  theta_0_2_1.c
//  
//
//  Created by Harry Jiang on 2016-08-16.
//
//

#include <math.h>
#include "mex.h"

mxArray* call_eval_fd(mxArray* rng, mxArray* fdobj, long Lfdobj){
    mxArray *result[1];
    int nin;
    
    if (Lfdobj == NULL) {
        mxArray *input[2] = {rng, fdobj};
        nin = 2;
        mexCallMATLAB(1, result, nin, input, "eval_fd");
    }
    else {
        mxArray *LFDOBJ     = mxCreateDoubleMatrix(1, 1, mxREAL);
        double *lfd         = mxGetPr(LFDOBJ);
        lfd[0]              = Lfdobj;
        mxArray *input[3]   = {rng, fdobj, LFDOBJ};
        
        nin = 3;
        mexCallMATLAB(1, result, nin, input, "eval_fd");
    }
    
    return result[0];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define HANEW_OUT       plhs[0]
    #define FAILINDEX_OUT   plhs[1]
    #define FAC_OUT         plhs[2]
    
    // Macros for input
    #define HANEW_IN        prhs[0]
    #define FAILINDEX_IN    prhs[1]
    #define FAC_IN          prhs[2]
    #define WFDOBJ          prhs[3]
    #define TINANEW         prhs[4]
    #define UA              prhs[5]
    #define SLAMBDA         prhs[6]
    #define TQIN            prhs[7]
    #define N_IN            prhs[8]
    #define HOLD            prhs[9]
    #define CRIT            prhs[10]
    #define ACTIVE          prhs[11]
    
    // Compute function values
    mxArray *WMATA  = call_eval_fd(TINANEW, WFDOBJ, NULL);
    double *wmata   = mxGetPr(WMATA);
    
    long M          = mxGetM(TINANEW);
    long N          = mxGetN(WMATA);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    double emata;
    
    HANEW_OUT       = HANEW_IN;
    double *hanew   = mxGetPr(HANEW_OUT);
    
    double *ua      = mxGetPr(UA);
    
    FAILINDEX_OUT   = FAILINDEX_IN;
    mxLogical *fail = mxGetLogicals(FAILINDEX_OUT);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    double *tinanew = mxGetPr(TINANEW);
    double *tqin    = mxGetPr(TQIN);
    double *active  = mxGetPr(ACTIVE);
    double *hold    = mxGetPr(HOLD);
    
    double crit     = mxGetPr(CRIT)[0];
    
    // Halve fac in preparation for next step size iteration
    FAC_OUT         = FAC_IN;
    double *fac     = mxGetPr(FAC_OUT);
    fac[0]         /= 2;
    
    long m, n;
    for (m = 0; m < M; m++){
        if (fail[m]){
            
            hanew[m] = 0.0;
            
            for (n = 0; n < N; n++){
                long index  = m + M*n;
                
                emata       = exp(wmata[index]);
                // Current vector of function values
                hanew[m]   -= (ua[index] * wmata[index] - log(1 + emata));
            }
            
            hanew[m] /= Nsubjects;
            
            if (slambda > 0) {
                hanew[m]   += slambda * pow(tinanew[m] - tqin[(long)active[m] - 1], 2) / Nsubjects;
            }
        }
        
        fail[m] = hold[(long)active[m] - 1] - hanew[m] < -crit;
    }
    
    return;
}