//
//  theta_0_1_0.c
//  
//
//  Created by Harry Jiang on 2016-08-10.
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
    #define STEP        plhs[1]
    //#define TEST_OUT    plhs[2]
    
    // Macros for input
    #define ITER_IN     prhs[0]
    #define TINA        prhs[1]
    #define WFDOBJ      prhs[2]
    #define N_IN        prhs[3]
    #define UA          prhs[4]
    #define SLAMBDA     prhs[5]
    #define TQIN        prhs[6]
    #define ACTIVE      prhs[7]
    
    ITER_OUT        = ITER_IN;
    double *iter    = mxGetPr(ITER_OUT);
    iter[0]        += 1;
    
    mxArray *WMATA      = call_eval_fd(TINA, WFDOBJ, NULL);
    //mexPrintf("wmata eval called ");
    mxArray *DWMATA     = call_eval_fd(TINA, WFDOBJ, 1);
    //mexPrintf("dwmata eval called ");
    mxArray *D2WMATA    = call_eval_fd(TINA, WFDOBJ, 2);
    //mexPrintf("d2wmata eval called\n");
    
    double *wmata   = mxGetPr(WMATA);
    double *dwmata  = mxGetPr(DWMATA);
    double *d2wmata = mxGetPr(D2WMATA);
    
    long M          = mxGetM(WMATA);
    long N          = mxGetN(WMATA);
    long Nsubjects  = (long)mxGetPr(N_IN)[0];
    
    //mexPrintf("M = %d\n", M);
    //mexPrintf("N = %d\n", N);
    //mexPrintf("Nsubjects = %d\n", Nsubjects);
    
    double emata, pmata, rmata;
    
    double *ua      = mxGetPr(UA);
    
    double dha, d2ha;
    
    STEP            = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *step    = mxGetPr(STEP);
    
    //mxArray *TQACT  = mxCreateDoubleMatrix(M, 1, mxREAL);
    //double *tqact   = mxGetPr(TQACT);
    
    //mxArray *DINC   = mxCreateDoubleMatrix(M, 1, mxREAL);
    //double *dinc    = mxGetPr(DINC);
    
    double *tina    = mxGetPr(TINA);
    double *tqin    = mxGetPr(TQIN);
    double *active  = mxGetPr(ACTIVE);
    
    double slambda  = mxGetPr(SLAMBDA)[0];
    
    long m, n;
    for (m = 0; m < M; m++) {
        
        dha     = 0.0;
        d2ha    = 0.0;
        
        //dinc[m] = 0.0;
        
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
            
            //dinc[m] = ((2 * slambda * (tina[m] - tqin[(long)active[m] - 1])) / Nsubjects );
        }
        //tqact[m] = tqin[(long)active[m] - 1];
        //mexPrintf("%f ", active[m]);
        //mexPrintf("%f ", d2ha[m]);
        
        // Find a step size that reduces all fuction values in the active set
        // Computer the negative of initial step for each theta
        step[m] = dha / d2ha;
        if (step[m] < -2) step[m] = -2;
        if (step[m] >  2) step[m] = 2;
        //mexPrintf("\n");
    }
    
    //TEST_OUT = DINC;
    
    return;
}