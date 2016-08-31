//
//  theta_0_3_0.c
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
    #define TINA        plhs[0]
    #define TIN_OUT     plhs[1]
    #define ACTIVE_OUT  plhs[2]
    #define HOLD_OUT    plhs[3]
    #define UA_OUT      plhs[4]
    //#define TEST_OUT    plhs[5]
    
    // Macros for input
    #define TINANEW     prhs[0]
    #define TIN_IN      prhs[1]
    #define ACTIVE_IN   prhs[2]
    #define HOLD_IN     prhs[3]
    #define UA_IN       prhs[4]
    #define UIN         prhs[5]
    #define WFDOBJ      prhs[6]
    #define SLAMBDA     prhs[7]
    #define HANEW       prhs[8]
    #define N_IN        prhs[9]
    #define TQIN        prhs[10]
    #define CRIT        prhs[11]
    
    mxArray *WMATA      = call_eval_fd(TINANEW, WFDOBJ, NULL);
    mxArray *DWMATA     = call_eval_fd(TINANEW, WFDOBJ, 1);
    mxArray *D2WMATA    = call_eval_fd(TINANEW, WFDOBJ, 2);
    double *wmata   = mxGetPr(WMATA);
    double *dwmata  = mxGetPr(DWMATA);
    double *d2wmata = mxGetPr(D2WMATA);
    
    double emata, pmata, rmata;
    
    double *ua      = mxGetPr(UA_IN);
    
    long M          = mxGetM(WMATA);
    long N          = mxGetN(WMATA);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    TIN_OUT         = TIN_IN;
    double *tin     = mxGetPr(TIN_OUT);
    double *tinanew = mxGetPr(TINANEW);
    double *tqin    = mxGetPr(TQIN);
    double *active  = mxGetPr(ACTIVE_IN);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    double *hanew   = mxGetPr(HANEW);
    
    HOLD_OUT        = HOLD_IN;
    double *hold    = mxGetPr(HOLD_OUT);
    double dha, d2ha;
    
    double crit     = mxGetPr(CRIT)[0];
    
    long *activenew = (long *)calloc(M, sizeof(long));
    long nactivenew = 0;
    
    // Either all function values reduced or 10 halvings completed
    // Save current function values
    
    long m, n;
    for (m = 0; m < M; m++) {
        
        tin[(long)active[m] - 1]    = tinanew[m];
        
        dha     = 0.0;
        d2ha    = 0.0;
        
        for (n = 0; n < N; n++) {
            long index  = m + M*n;
            
            emata       = exp(wmata[index]);
            pmata       = emata/(1 + emata);
            rmata       = ua[index] - pmata;
            
            dha        += rmata * dwmata[index];
            d2ha       += ( (-d2wmata[index] * rmata) +
                           (pow(dwmata[index], 2) * pmata * (1 - pmata)) );
        }
        
        dha     = -dha/Nsubjects;
        d2ha    = d2ha/Nsubjects;
        
        if (slambda > 0) {
            dha    += (2 * slambda * (tinanew[m] - tqin[(long)active[m] - 1])) / Nsubjects;
            d2ha   += 2 * slambda / Nsubjects;
        }
        //mexPrintf("%f ", dha);
        //mexPrintf("%f\n", d2ha);
        
        hold[(long)active[m] - 1]   = hanew[m];
        
        if (fabs(dha) > crit && d2ha > 0) {
            //mexPrintf("%d\n", nactivenew);
            activenew[nactivenew++] = (long)active[m];
        }
    }
    //mexPrintf("%d\n", nactivenew);
    
    ACTIVE_OUT          = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *activeout   = mxGetPr(ACTIVE_OUT);
    
    TINA                = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *tina        = mxGetPr(TINA);
    
    UA_OUT              = mxCreateDoubleMatrix(nactivenew, N, mxREAL);
    double *ua_out      = mxGetPr(UA_OUT);
    double *uin         = mxGetPr(UIN);
    
    for (m = 0; m < nactivenew; m++) {
        
        activeout[m]    = activenew[m];
        tina[m]         = tin[activenew[m] - 1];
        
        //mexPrintf("%d    ", activenew[m]);
        
        for (n = 0; n < N; n++) {
            long newindex       = m + nactivenew*n;
            long index          = activenew[m] - 1 + Nsubjects*n;
            
            //mexPrintf("%1.0f    ", uin[index]);
            ua_out[newindex]    = uin[index];
        }
        //mexPrintf("\n");
    }
    
    free(activenew);
    //TEST_OUT = UIN;
    
    // Return for a new optimization step
    
    return;
}