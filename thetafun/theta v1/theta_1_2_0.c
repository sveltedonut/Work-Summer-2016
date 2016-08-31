//
//  theta_1_2_0.c
//  
//
//  Created by Harry Jiang on 2016-08-18.
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
    #define STEPITER_OUT    plhs[0]
    #define TINANEW_OUT     plhs[1]
    #define HANEW_OUT       plhs[2]
    #define FAILINDEX_OUT   plhs[3]
    #define FAC_OUT         plhs[4]
    
    // Macros for input
    #define STEPITER_IN     prhs[0]
    #define STEP            prhs[1]
    #define SCRMIN          prhs[2]
    #define SCRMAX          prhs[3]
    #define FAC             prhs[4]
    #define TINANEW_IN      prhs[5]
    #define ACTIVE          prhs[6]
    #define HANEW_IN        prhs[7]
    #define FAILINDEX_IN    prhs[8]
    #define FAC_IN          prhs[9]
    #define WFDOBJ          prhs[10]
    #define UA              prhs[11]
    #define SLAMBDA         prhs[12]
    #define TQIN            prhs[13]
    #define N_IN            prhs[14]
    #define HOLD            prhs[15]
    #define CRIT            prhs[16]
    
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
    
    FAC_OUT         = FAC_IN;
    double *fac     = mxGetPr(FAC_OUT);
    
    FAILINDEX_OUT   = FAILINDEX_IN;
    mxLogical *fail = mxGetLogicals(FAILINDEX_OUT);
    
    long f = 0;
    long m;
    for (m = 0; m < M; m++) {
        // Compute current step size
        if (fail[m]) {
            stepnew     = fac[0] * step[m];
            tinanewfail = tinanew[m];
            if (tinanewfail - stepnew <= scrmin) stepnew = tinanewfail - scrmin * 1.01;
            if (tinanewfail - stepnew >= scrmax) stepnew = -(scrmax * 0.99 - tinanewfail);
            
            // Update ability values
            tinanew[m] -= stepnew;
            
            f++;
        }
    }
    
    // Compute function values
    mxArray *WMATA  = call_eval_fd(TINANEW_OUT, WFDOBJ, NULL);
    double *wmata   = mxGetPr(WMATA);
    
    long N          = mxGetN(WMATA);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    double emata;
    
    HANEW_OUT       = HANEW_IN;
    double *hanew   = mxGetPr(HANEW_OUT);
    
    double *ua      = mxGetPr(UA);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    double *tqin    = mxGetPr(TQIN);
    double *hold    = mxGetPr(HOLD);
    
    double crit     = mxGetPr(CRIT)[0];
    
    // Halve fac in preparation for next step size iteration
    fac[0]         /= 2;
    
    long n;
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