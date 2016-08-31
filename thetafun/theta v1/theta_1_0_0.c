//
//  theta_1_0_0.c
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
    #define ITER        plhs[0]
    #define TINA        plhs[1]
    #define ACTIVE      plhs[2]
    #define UA          plhs[3]
    #define HOLD        plhs[4]
    
    // Macros for input
    #define TIN         prhs[0]
    #define WFDOBJ      prhs[1]
    #define UIN         prhs[2]
    #define N_IN        prhs[3]
    #define SLAMBDA     prhs[4]
    #define TQIN        prhs[5]
    #define CRIT        prhs[6]
    
    mxArray *WMAT   = call_eval_fd(TIN, WFDOBJ, NULL);
    mxArray *DWMAT  = call_eval_fd(TIN, WFDOBJ, 1);
    mxArray *D2WMAT = call_eval_fd(TIN, WFDOBJ, 2);
    double *wmat    = mxGetPr(WMAT);
    double *dwmat   = mxGetPr(DWMAT);
    double *d2wmat  = mxGetPr(D2WMAT);
    
    double emat, pmat, rmat;
    
    double *uin     = mxGetPr(UIN);
    
    long M          = mxGetM(WMAT);
    long N          = mxGetN(WMAT);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    double *tin     = mxGetPr(TIN);
    double *tqin    = mxGetPr(TQIN);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    HOLD            = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *hold    = mxGetPr(HOLD);
    double dha, d2ha;
    
    double crit     = mxGetPr(CRIT)[0];
    
    long *activenew = (long *)calloc(M, sizeof(long));
    long nactivenew = 0;
    
    // Initialize iterations
    
    ITER            = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *iter    = mxGetPr(ITER);
    iter[0]         = 0;
    
    long m, n;
    for (m = 0; m < M; m++) {
        
        hold[m] = 0.0;
        dha     = 0.0;
        d2ha    = 0.0;
        
        for (n = 0; n < N; n++) {
            long index  = m + M*n;
            
            // Initial function, first and second derivative values  wrt theta
            
            emat        = exp(wmat[index]);
            pmat        = emat/(1 + emat);
            rmat        = uin[index] - pmat;
            
            hold[m]    += uin[index] * wmat[index] - log(1 + emat);
            dha        += rmat * dwmat[index];
            d2ha       += ( (-d2wmat[index] * rmat) +
                           (pow(dwmat[index], 2) * pmat * (1 - pmat)) );
        }
        hold[m] = -hold[m]/Nsubjects;
        dha     = -dha/Nsubjects;
        d2ha    = d2ha/Nsubjects;
        
        if (slambda > 0) {
            hold[m]    += slambda * pow(tin[m] - tqin[m], 2) / Nsubjects;
            dha        += (2 * slambda * (tin[m] - tqin[m])) / Nsubjects;
            d2ha       += 2 * slambda / Nsubjects;
        }
        
        // Find cases where gradient size exceeds criterion and further optimization is required
        
        if (fabs(dha) > crit && d2ha > 0) {
            activenew[nactivenew++] = m + 1;
        }
    }
    
    ACTIVE          = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *active  = mxGetPr(ACTIVE);
    
    TINA            = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *tina    = mxGetPr(TINA);
    
    UA              = mxCreateDoubleMatrix(nactivenew, N, mxREAL);
    double *ua      = mxGetPr(UA);
    
    for (m = 0; m < nactivenew; m++) {
        
        active[m]   = activenew[m];
        
        // Select active theta's and data
        
        tina[m]     = tin[activenew[m] - 1];
        
        //mexPrintf("%d    ", activenew[m]);
        
        for (n = 0; n < N; n++) {
            long newindex   = m + nactivenew*n;
            long index      = activenew[m] - 1 + Nsubjects*n;
            
            //mexPrintf("%1.0f    ", uin[index]);
            ua[newindex]    = uin[index];
        }
        //mexPrintf("\n");
    }
    
    free(activenew);
    
    // Update scores in active set by a Newton-Raphson step
    
    return;
}