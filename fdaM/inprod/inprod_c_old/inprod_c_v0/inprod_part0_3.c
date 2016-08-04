//
//  inprod_part0_3.c
//  
//
//  Created by Harry Jiang on 2016-05-10.
//
//

#include "mex.h"

// call isa_fd from matlab
bool call_isa_fd(mxArray* fdobj) {
    mxArray *result[1];
    bool *isa_fd;
    
    mexCallMATLAB(1, result, 1, &fdobj, "isa_fd");
    
    isa_fd = (bool*)mxGetPr(result[0]);
    //mexPrintf("is a fd: %d\n", isa_fd[0]);
    
    return isa_fd[0];
}

// call eval_fd from matlab
mxArray* call_eval_fd(mxArray* rng, mxArray* fdobj, mxArray* Lfdobj){
    mxArray *result[1];
    int nin;
    
    if (Lfdobj == NULL) {
        mxArray *input[2] = {rng, fdobj};
        nin = 2;
        mexCallMATLAB(1, result, nin, input, "eval_fd");
    }
    else {
        mxArray *input[3] = {rng, fdobj, Lfdobj};
        nin = 3;
        mexCallMATLAB(1, result, nin, input, "eval_fd");
    }
    
    return result[0];
}

//call isempty from matlab
bool call_isempty(mxArray* wtfd) {
    mxArray *result[1];
    bool *isempty;
    
    mexCallMATLAB(1, result, 1, &wtfd, "isempty");
    
    isempty = (bool*)mxGetPr(result[0]);
    
    return isempty[0];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define FX1         plhs[0]
    #define FX2         plhs[1]
    
    // Macros for input
    #define RNGI        prhs[0]
    #define FDOBJ1      prhs[1]
    #define FDOBJ2      prhs[2]
    #define LFDOBJ1     prhs[3]
    #define LFDOBJ2     prhs[4]
    #define WTFD        prhs[5]
    #define NREP2       prhs[6]
    
    
    if (call_isa_fd(FDOBJ1)) {
        FX1 = call_eval_fd(RNGI, FDOBJ1, LFDOBJ1);
    }
    
    if (call_isa_fd(FDOBJ2)) {
        FX2 = call_eval_fd(RNGI, FDOBJ2, LFDOBJ2);
    }
    
    //apply weight to fx2
    if (call_isempty(WTFD) == false) {
        mxArray *WTD    = call_eval_fd(RNGI, WTFD, NULL);
        double *wtd     = mxGetPr(WTD);
        double *fx2     = mxGetPr(FX2);
        int M           = mxGetM(WTD);
        int N           = mxGetN(WTD);
        int i, j;
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                fx2[M*i + j] *= wtd[i];
            }
        }
    }
    
    return;
}