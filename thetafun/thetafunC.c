//
//  theta_6_0_0.c
//  
//
//  Created by Harry Jiang on 2016-08-18.
//
//

#include <math.h>
#include "mex.h"
#include "bsplineC.h"

// Main function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define ITER        plhs[0]
    #define TIN_OUT     plhs[1]
    #define HOLD        plhs[2]
    #define DHOLD       plhs[3]
    #define D2HOLD      plhs[4]
    
    // Macros for input
    #define ITERMAX     prhs[0]
    #define TIN_IN      prhs[1]
    #define WBREAK      prhs[2]
    #define WORDER      prhs[3]
    #define WCOEFS      prhs[4]
    #define UIN         prhs[5]
    #define N_IN        prhs[6]
    #define SLAMBDA     prhs[7]
    #define TQIN        prhs[8]
    #define SCRRNG      prhs[9]
    #define CRIT        prhs[10]
    
    mxArray *WMAT   = call_eval_fd(TIN_IN, WBREAK, WORDER, WCOEFS, 0);
    mxArray *DWMAT  = call_eval_fd(TIN_IN, WBREAK, WORDER, WCOEFS, 1);
    mxArray *D2WMAT = call_eval_fd(TIN_IN, WBREAK, WORDER, WCOEFS, 2);
    double *wmat    = mxGetPr(WMAT);
    double *dwmat   = mxGetPr(DWMAT);
    double *d2wmat  = mxGetPr(D2WMAT);
    
    double emat, pmat, rmat;
    
    double *uin     = mxGetPr(UIN);
    
    long M          = mxGetM(WMAT);
    long N          = mxGetN(WMAT);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    double *tin     = mxGetPr(TIN_IN);
    double *tqin    = mxGetPr(TQIN);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    HOLD            = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *hold    = mxGetPr(HOLD);
    DHOLD           = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *dhold   = mxGetPr(DHOLD);
    D2HOLD          = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *d2hold  = mxGetPr(D2HOLD);
    
    double crit     = mxGetPr(CRIT)[0];
    
    long *activenew = (long *)calloc(M, sizeof(long));
    long nactivenew = 0;
    
    // Compute initial probability values for interior scores
    long m, n;
    for (m = 0; m < M; m++) {
        
        // Initial function, first and second derivative values  wrt theta
        hold[m]     = 0.0;
        dhold[m]    = 0.0;
        d2hold[m]   = 0.0;
        
        for (n = 0; n < N; n++) {
            long index  = m + M*n;
            
            emat        = exp(wmat[index]);
            pmat        = emat/(1 + emat);
            rmat        = uin[index] - pmat;
            
            hold[m]    += uin[index] * wmat[index] - log(1 + emat);
            dhold[m]   += rmat * dwmat[index];
            d2hold[m]  += ( (-d2wmat[index] * rmat) +
                           (pow(dwmat[index], 2) * pmat * (1 - pmat)) );
        }
        hold[m]     = -hold[m]/Nsubjects;
        dhold[m]    = -dhold[m]/Nsubjects;
        d2hold[m]   = d2hold[m]/Nsubjects;
        
        if (slambda > 0) {
            hold[m]    += slambda * pow(tin[m] - tqin[m], 2) / Nsubjects;
            dhold[m]   += (2 * slambda * (tin[m] - tqin[m])) / Nsubjects;
            d2hold[m]  += 2 * slambda / Nsubjects;
        }
        
        // Find cases where gradient size exceeds criterion and further optimization is required
        
        if (fabs(dhold[m]) > crit && d2hold[m] > 0) {
            activenew[nactivenew++] = m + 1;
        }
    }
    
    mxArray *ACTIVE = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *active  = mxGetPr(ACTIVE);
    
    mxArray *TINA   = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
    double *tina    = mxGetPr(TINA);
    
    mxArray *UA     = mxCreateDoubleMatrix(nactivenew, N, mxREAL);
    double *ua      = mxGetPr(UA);
    
    for (m = 0; m < nactivenew; m++) {
        
        active[m]   = activenew[m];
        
        // Select active theta's and data
        
        tina[m]     = tin[activenew[m] - 1];
        
        for (n = 0; n < N; n++) {
            long newindex   = m + nactivenew*n;
            long index      = activenew[m] - 1 + Nsubjects*n;
            
            ua[newindex]    = uin[index];
        }
    }
    
    free(activenew);
    
    // Update scores in active set by a Newton-Raphson step
    
    ITER            = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *iter    = mxGetPr(ITER);
    iter[0]         = 0;
    long itermax    = (long)mxGetPr(ITERMAX)[0];
    
    nactivenew = 1;
    
    TIN_OUT         = TIN_IN;
    
    while (nactivenew > 0 && iter[0] < itermax) {
        iter[0]++;
        
        WMAT            = call_eval_fd(TINA, WBREAK, WORDER, WCOEFS, 0);
        DWMAT           = call_eval_fd(TINA, WBREAK, WORDER, WCOEFS, 1);
        D2WMAT          = call_eval_fd(TINA, WBREAK, WORDER, WCOEFS, 2);
        wmat            = mxGetPr(WMAT);
        dwmat           = mxGetPr(DWMAT);
        d2wmat          = mxGetPr(D2WMAT);
        
        M               = mxGetM(WMAT);
        N               = mxGetN(WMAT);
        
        ua              = mxGetPr(UA);
        
        double *step    = (double *)malloc(M * sizeof(double));
        double stepnew;
        
        tina            = mxGetPr(TINA);
        tqin            = mxGetPr(TQIN);
        active          = mxGetPr(ACTIVE);
        
        double scrmin   = mxGetPr(SCRRNG)[0];
        double scrmax   = mxGetPr(SCRRNG)[1];
        
        tin             = mxGetPr(TIN_IN);
        
        mxArray *TINANEW    = mxCreateDoubleMatrix(M, 1, mxREAL);
        double *tinanew     = mxGetPr(TINANEW);
        double tinanewfail;
        
        double *hanew   = (double *)malloc(M * sizeof(double));
        
        // Initial set of thetas failing to meet criterion
        mxArray *FAILINDEX  = mxCreateLogicalMatrix(M, 1);
        mxLogical *fail     = mxGetLogicals(FAILINDEX);
        
        long m, n;
        for (m = 0; m < M; m++) {
            long activeindex = (long)active[m] - 1;
            
            dhold[activeindex]    = 0.0;
            d2hold[activeindex]   = 0.0;
            
            // Compute probabilites and other quantities for active cases
            for (n = 0; n < N; n++) {
                long index  = m + M*n;
                
                emat       = exp(wmat[index]);
                pmat       = emat/(1 + emat);
                rmat       = ua[index] - pmat;
                
                dhold[activeindex]   += rmat * dwmat[index];
                d2hold[activeindex]  += ( (-d2wmat[index] * rmat) +
                                         (pow(dwmat[index], 2) * pmat * (1 - pmat)) );
            }
            
            dhold[activeindex]    = -dhold[activeindex]/Nsubjects;
            d2hold[activeindex]   = d2hold[activeindex]/Nsubjects;
            
            // Add regularization terms if needed
            if (slambda > 0) {
                dhold[activeindex]  += (2 * slambda * (tina[m] - tqin[activeindex])) / Nsubjects;
                d2hold[activeindex] += 2 * slambda / Nsubjects;
            }
            
            // Find a step size that reduces all fuction values in the active set
            // Computer the negative of initial step for each theta
            step[m] = dhold[activeindex] / d2hold[activeindex];
            if (step[m] < -2) step[m] = -2;
            if (step[m] >  2) step[m] = 2;
            
            // Ensure that initial step does not go to boundaries
            if (tina[m] - step[m] <= scrmin) step[m] = tina[m] - scrmin * 1.01;
            if (tina[m] - step[m] >= scrmax) step[m] = -(scrmax * 0.99 - tina[m]);
            
            tinanew[m]  = tin[activeindex];
            hanew[m]    = hold[activeindex];
            fail[m]     = true;
        }
        
        // Half fac as required to achieve the reduction of all function values
        double fac      = 1;
        
        long loopfail   = 1;
        
        long stepiter   = 0;
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
            WMAT    = call_eval_fd(TINANEW, WBREAK, WORDER, WCOEFS, 0);
            wmat    = mxGetPr(WMAT);
            N       = mxGetN(WMAT);
            
            // Halve fac in preparation for next step size iteration
            fac    /= 2;
            
            long n;
            for (m = 0; m < M; m++){
                long activeindex = (long)active[m] - 1;
                
                if (fail[m]){
                    
                    hanew[m] = 0.0;
                    
                    for (n = 0; n < N; n++){
                        long index  = m + M*n;
                        
                        emat       = exp(wmat[index]);
                        // Current vector of function values
                        hanew[m]   -= (ua[index] * wmat[index] - log(1 + emat));
                    }
                    
                    hanew[m] /= Nsubjects;
                    
                    if (slambda > 0) {
                        hanew[m]   += slambda * pow(tinanew[m] - tqin[activeindex], 2) / Nsubjects;
                    }
                }
                
                fail[m] = hold[activeindex] - hanew[m] < -crit;
                
                if (fail[m]) loopfail++;
            }
        }
        
        WMAT           = call_eval_fd(TINANEW, WBREAK, WORDER, WCOEFS, 0);
        DWMAT          = call_eval_fd(TINANEW, WBREAK, WORDER, WCOEFS, 1);
        D2WMAT         = call_eval_fd(TINANEW, WBREAK, WORDER, WCOEFS, 2);
        wmat   = mxGetPr(WMAT);
        dwmat  = mxGetPr(DWMAT);
        d2wmat = mxGetPr(D2WMAT);
        
        long *activenew = (long *)calloc(M, sizeof(long));
        nactivenew      = 0;
        
        // Either all function values reduced or 10 halvings completed
        // Save current function values
        
        for (m = 0; m < M; m++) {
            long activeindex = (long)active[m] - 1;
            
            tin[activeindex]    = tinanew[m];
            
            dhold[activeindex]    = 0.0;
            d2hold[activeindex]   = 0.0;
            
            for (n = 0; n < N; n++) {
                long index  = m + M*n;
                
                emat       = exp(wmat[index]);
                pmat       = emat/(1 + emat);
                rmat       = ua[index] - pmat;
                
                dhold[activeindex]   += rmat * dwmat[index];
                d2hold[activeindex]  += ( (-d2wmat[index] * rmat) +
                                         (pow(dwmat[index], 2) * pmat * (1 - pmat)) );
            }
            
            dhold[activeindex]    = -dhold[activeindex]/Nsubjects;
            d2hold[activeindex]   = d2hold[activeindex]/Nsubjects;
            
            if (slambda > 0) {
                dhold[activeindex]  += (2 * slambda * (tinanew[m] - tqin[activeindex])) / Nsubjects;
                d2hold[activeindex] += 2 * slambda / Nsubjects;
            }
            
            hold[activeindex]   = hanew[m];
            
            if (fabs(dhold[activeindex]) > crit && d2hold[activeindex] > 0) {
                activenew[nactivenew++] = (long)active[m];
            }
        }
        ACTIVE              = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
        double *activeout   = mxGetPr(ACTIVE);
        
        TINA                = mxCreateDoubleMatrix(nactivenew, 1, mxREAL);
        double *tina_out    = mxGetPr(TINA);
        
        UA                  = mxCreateDoubleMatrix(nactivenew, N, mxREAL);
        double *ua_out      = mxGetPr(UA);
        double *uin         = mxGetPr(UIN);
        
        for (m = 0; m < nactivenew; m++) {
            
            activeout[m]    = activenew[m];
            tina_out[m]     = tin[activenew[m] - 1];
            
            for (n = 0; n < N; n++) {
                long newindex       = m + nactivenew*n;
                long index          = activenew[m] - 1 + Nsubjects*n;
                
                ua_out[newindex]    = uin[index];
            }
        }
        
        free(hanew);
        free(step);
        free(activenew);
        
    }
    
    return;
    
}
