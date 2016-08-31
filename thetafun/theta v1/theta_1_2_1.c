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
    #define TINANEW_OUT     plhs[0]
    #define HANEW_OUT       plhs[1]
    
    // Macros for input
    #define STEP            prhs[0]
    #define SCRMIN          prhs[1]
    #define SCRMAX          prhs[2]
    #define TINANEW_IN      prhs[3]
    #define ACTIVE          prhs[4]
    #define HANEW_IN        prhs[5]
    #define FAILINDEX       prhs[6]
    #define WFDOBJ          prhs[7]
    #define UA              prhs[8]
    #define SLAMBDA         prhs[9]
    #define TQIN            prhs[10]
    #define N_IN            prhs[11]
    #define HOLD            prhs[12]
    #define CRIT            prhs[13]
    
    long M              = mxGetM(STEP);
    
    double *step        = mxGetPr(STEP);
    double stepnew;
    
    double *active      = mxGetPr(ACTIVE);
    
    TINANEW_OUT         = TINANEW_IN;
    double *tinanew     = mxGetPr(TINANEW_OUT);
    double tinanewfail;
    
    double scrmin       = mxGetPr(SCRMIN)[0];
    double scrmax       = mxGetPr(SCRMAX)[0];
    
    // Half fac as required to achieve the reduction of all function values
    double fac          = 1;
    
    mxLogical *fail     = mxGetLogicals(FAILINDEX);
    long loopfail = 1;
    
    long Nsubjects      = (long)mxGetPr(N_IN)[0];
    
    double emata;
    
    HANEW_OUT           = HANEW_IN;
    double *hanew       = mxGetPr(HANEW_OUT);
    
    double *ua          = mxGetPr(UA);
    
    double slambda      = mxGetPr(SLAMBDA)[0];
    
    double *tqin        = mxGetPr(TQIN);
    double *hold        = mxGetPr(HOLD);
    
    double crit         = mxGetPr(CRIT)[0];
    
    long stepiter       = 0;
    while (loopfail > 0 && stepiter < 10){
        
        stepiter++;
        loopfail    = 0;
        
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
            }
        }
        
        // Compute function values
        mxArray *WMATA  = call_eval_fd(TINANEW_OUT, WFDOBJ, NULL);
        double *wmata   = mxGetPr(WMATA);
        
        long N          = mxGetN(WMATA);
        
        // Halve fac in preparation for next step size iteration
        fac            /= 2;
        
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
            
            if (fail[m]) loopfail++;        }
        
    }
    
    return;
}