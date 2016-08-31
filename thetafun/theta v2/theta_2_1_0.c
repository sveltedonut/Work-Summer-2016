//
//  theta_2_1_0.c
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
    #define TINANEW     plhs[1]
    #define HANEW       plhs[2]
    
    // Macros for input
    #define ITER_IN     prhs[0]
    #define TINA        prhs[1]
    #define WFDOBJ      prhs[2]
    #define N_IN        prhs[3]
    #define UA          prhs[4]
    #define SLAMBDA     prhs[5]
    #define TQIN        prhs[6]
    #define ACTIVE      prhs[7]
    #define SCRMIN      prhs[8]
    #define SCRMAX      prhs[9]
    #define TIN         prhs[10]
    #define HOLD        prhs[11]
    #define CRIT        prhs[12]
    
    ITER_OUT        = ITER_IN;
    double *iter    = mxGetPr(ITER_OUT);
    iter[0]        += 1;
    
    mxArray *WMATA      = call_eval_fd(TINA, WFDOBJ, NULL);
    mxArray *DWMATA     = call_eval_fd(TINA, WFDOBJ, 1);
    mxArray *D2WMATA    = call_eval_fd(TINA, WFDOBJ, 2);
    
    double *wmata   = mxGetPr(WMATA);
    double *dwmata  = mxGetPr(DWMATA);
    double *d2wmata = mxGetPr(D2WMATA);
    
    long M          = mxGetM(WMATA);
    long N          = mxGetN(WMATA);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    double emata, pmata, rmata;
    
    double *ua      = mxGetPr(UA);
    
    double dha, d2ha;
    
    double *step    = (double *)malloc(M * sizeof(double));
    double stepnew;
    
    double *tina    = mxGetPr(TINA);
    double *tqin    = mxGetPr(TQIN);
    double *active  = mxGetPr(ACTIVE);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    double scrmin   = mxGetPr(SCRMIN)[0];
    double scrmax   = mxGetPr(SCRMAX)[0];
    
    double *tin     = mxGetPr(TIN);
    double *hold    = mxGetPr(HOLD);
    
    TINANEW         = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *tinanew = mxGetPr(TINANEW);
    double tinanewfail;
    
    HANEW           = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *hanew   = mxGetPr(HANEW);
    
    // Initial set of thetas failing to meet criterion
    mxArray *FAILINDEX  = mxCreateLogicalMatrix(M, 1);
    mxLogical *fail     = mxGetLogicals(FAILINDEX);
    
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
    
    // Half fac as required to achieve the reduction of all function values
    double fac          = 1;
    
    long loopfail = 1;
    
    double crit         = mxGetPr(CRIT)[0];
    
    long stepiter       = 0;while (loopfail > 0 && stepiter < 10){
        
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
        mxArray *WMATA  = call_eval_fd(TINANEW, WFDOBJ, NULL);
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