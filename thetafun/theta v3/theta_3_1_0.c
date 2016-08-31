//
//  theta_3_1_0.c
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
    #define ITER_OUT    plhs[0]
    #define TINA_OUT    plhs[1]
    #define TIN_OUT     plhs[2]
    #define ACTIVE_OUT  plhs[3]
    #define HOLD_OUT    plhs[4]
    #define UA_OUT      plhs[5]
    
    // Macros for input
    #define ITER_IN     prhs[0]
    #define TINA_IN     prhs[1]
    #define TIN_IN      prhs[2]
    #define ACTIVE_IN   prhs[3]
    #define HOLD_IN     prhs[4]
    #define UA_IN       prhs[5]
    #define UIN         prhs[6]
    #define WFDOBJ      prhs[7]
    #define N_IN        prhs[8]
    #define SLAMBDA     prhs[9]
    #define TQIN        prhs[10]
    #define SCRMIN      prhs[11]
    #define SCRMAX      prhs[12]
    #define CRIT        prhs[13]

    ITER_OUT        = ITER_IN;
    double *iter    = mxGetPr(ITER_OUT);
    iter[0]        += 1;
    
    mxArray *WMATA      = call_eval_fd(TINA_IN, WFDOBJ, NULL);
    mxArray *DWMATA     = call_eval_fd(TINA_IN, WFDOBJ, 1);
    mxArray *D2WMATA    = call_eval_fd(TINA_IN, WFDOBJ, 2);
    
    double *wmata   = mxGetPr(WMATA);
    double *dwmata  = mxGetPr(DWMATA);
    double *d2wmata = mxGetPr(D2WMATA);
    
    long M          = mxGetM(WMATA);
    long N          = mxGetN(WMATA);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    double emata, pmata, rmata;
    
    double *ua      = mxGetPr(UA_IN);
    
    double *step    = (double *)malloc(M * sizeof(double));
    double stepnew;
    
    double *tina    = mxGetPr(TINA_IN);
    double *tqin    = mxGetPr(TQIN);
    double *active  = mxGetPr(ACTIVE_IN);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    double scrmin   = mxGetPr(SCRMIN)[0];
    double scrmax   = mxGetPr(SCRMAX)[0];
    
    TIN_OUT         = TIN_IN;
    double *tin     = mxGetPr(TIN_IN);
    
    HOLD_OUT        = HOLD_IN;
    double *hold    = mxGetPr(HOLD_OUT);
    double dha, d2ha;
    
    mxArray *TINANEW    = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *tinanew     = mxGetPr(TINANEW);
    double tinanewfail;
    
    double *hanew   = (double *)malloc(M * sizeof(double));
    
    // Initial set of thetas failing to meet criterion
    mxArray *FAILINDEX  = mxCreateLogicalMatrix(M, 1);
    mxLogical *fail     = mxGetLogicals(FAILINDEX);
    
    //mexPrintf("part a\n");
    
    long m, n;
    for (m = 0; m < M; m++) {
        
        dha     = 0.0;
        d2ha    = 0.0;
        
        // Compute probabilites and other quantities for active cases
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
        
        // Add regularization terms if needed
        if (slambda > 0) {
            dha    += (2 * slambda * (tina[m] - tqin[(long)active[m] - 1])) / Nsubjects;
            d2ha   += 2 * slambda / Nsubjects;
        }
        
        // Find a step size that reduces all fuction values in the active set
        // Computer the negative of initial step for each theta
        step[m] = dha / d2ha;
        if (step[m] < -2) step[m] = -2;
        if (step[m] >  2) step[m] = 2;
        
        // Ensure that initial step does not go to boundaries
        if (tina[m] - step[m] <= scrmin) step[m] = tina[m] - scrmin * 1.01;
        if (tina[m] - step[m] >= scrmax) step[m] = -(scrmax * 0.99 - tina[m]);
        
        tinanew[m]  = tin[(long)active[m] - 1];
        hanew[m]    = hold[(long)active[m] - 1];
        fail[m]     = true;
    }
    //mexPrintf("part b\n");
    
    // Half fac as required to achieve the reduction of all function values
    double fac          = 1;
    
    long loopfail = 1;
    
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
        WMATA           = call_eval_fd(TINANEW, WFDOBJ, NULL);
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
    //mexPrintf("part c\n");
    
    WMATA           = call_eval_fd(TINANEW, WFDOBJ, NULL);
    DWMATA          = call_eval_fd(TINANEW, WFDOBJ, 1);
    D2WMATA         = call_eval_fd(TINANEW, WFDOBJ, 2);
    wmata   = mxGetPr(WMATA);
    dwmata  = mxGetPr(DWMATA);
    d2wmata = mxGetPr(D2WMATA);
    
    long *activenew = (long *)calloc(M, sizeof(long));
    long nactivenew = 0;
    //mexPrintf("part d\n");
    
    // Either all function values reduced or 10 halvings completed
    // Save current function values
    
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
    ACTIVE_OUT          = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *activeout   = mxGetPr(ACTIVE_OUT);
    
    TINA_OUT            = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *tina_out    = mxGetPr(TINA_OUT);
    
    UA_OUT              = mxCreateDoubleMatrix(nactivenew, N, mxREAL);
    double *ua_out      = mxGetPr(UA_OUT);
    double *uin         = mxGetPr(UIN);
    
    for (m = 0; m < nactivenew; m++) {
        
        activeout[m]    = activenew[m];
        tina_out[m]     = tin[activenew[m] - 1];
        
        //mexPrintf("%d    ", activenew[m]);
        
        for (n = 0; n < N; n++) {
            long newindex       = m + nactivenew*n;
            long index          = activenew[m] - 1 + Nsubjects*n;
            
            //mexPrintf("%1.0f    ", uin[index]);
            ua_out[newindex]    = uin[index];
        }
        //mexPrintf("\n");
    }
    //mexPrintf("part e\n");
    
    free(hanew);
    free(step);
    free(activenew);
    
    return;
    
}
