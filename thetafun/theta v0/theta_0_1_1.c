//
//  theta_0_1_1.c
//  
//
//  Created by Harry Jiang on 2016-08-15.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define TINANEW     plhs[0]
    #define HANEW       plhs[1]
    #define FAC         plhs[2]
    #define STEPITER    plhs[3]
    #define FAILINDEX   plhs[4]
    
    // Macros for input
    #define TINA_IN     prhs[0]
    #define STEP        prhs[1]
    #define SCRMIN      prhs[2]
    #define SCRMAX      prhs[3]
    #define T_IN        prhs[4]
    #define ACTIVE      prhs[5]
    #define HOLD        prhs[6]
    
    long M          = mxGetM(STEP);
    
    double *step    = mxGetPr(STEP);
    
    double *tina    = mxGetPr(TINA_IN);
    
    double scrmin   = mxGetPr(SCRMIN)[0];
    double scrmax   = mxGetPr(SCRMAX)[0];
    
    double *active  = mxGetPr(ACTIVE);
    double *tin     = mxGetPr(T_IN);
    double *hold    = mxGetPr(HOLD);
    
    TINANEW         = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *tinanew = mxGetPr(TINANEW);
    
    HANEW           = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *hanew   = mxGetPr(HANEW);
    
    // Half fac as required to achieve the reduction of all function values
    FAC             = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *fac     = mxGetPr(FAC);
    fac[0]          = 1;
    
    STEPITER        = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *siter   = mxGetPr(STEPITER);
    siter[0]        = 0;
    
    // Initial set of thetas failing to meet criterion
    FAILINDEX       = mxCreateLogicalMatrix(M, 1);
    mxLogical *fail = mxGetLogicals(FAILINDEX);
    
    long m;
    for (m = 0; m < M; m++) {
        // Ensure that initial step does not go to boundaries
        if (tina[m] - step[m] <= scrmin) step[m] = tina[m] - scrmin * 1.01;
        if (tina[m] - step[m] >= scrmax) step[m] = -(scrmax * 0.99 - tina[m]);
        
        tinanew[m]  = tin[(long)active[m] - 1];
        hanew[m]    = hold[(long)active[m] - 1];
        fail[m]     = true;
    }
    
    return;
}