//
//  inprod_part1_2.c
//  
//
//  Created by Harry Jiang on 2016-05-19.
//
//

#include "mex.h"

//call isa_fd from matlab
bool call_isa_fd(const mxArray* fdobj) {
    mxArray *result[1];
    bool *isa_fd;
    
    mexCallMATLAB(1, result, 1, &fdobj, "isa_fd");
    
    isa_fd = (bool*)mxGetPr(result[0]);
    //mexPrintf("is a fd: %d\n", isa_fd[0]);
    
    return isa_fd[0];
}

// call eval_fd from matlab
mxArray* call_eval_fd(mxArray* rng, const mxArray* fdobj, const mxArray* Lfdobj){
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

// call isempty from matlab
bool call_isempty(const mxArray* wtfd) {
    mxArray *result[1];
    bool *isempty;
    
    mexCallMATLAB(1, result, 1, &wtfd, "isempty");
    
    isempty = (bool*)mxGetPr(result[0]);
    
    return isempty[0];
}

// evaluate estimated solutions fx1 and fx2 of inner product
void eval_fds(mxArray** fx1,            mxArray** fx2,
              const mxArray*  fdobj1,   const mxArray* fdobj2,
              const mxArray*  Lfdobj1,  const mxArray* Lfdobj2,
              const mxArray*  wtfd,     mxArray* rngi) {
    
    // eval_fd for fx1 and fx2 and necessary checks
    if (call_isa_fd(fdobj1)) {
        *fx1 = call_eval_fd(rngi, fdobj1, Lfdobj1);
    }
    
    if (call_isa_fd(fdobj2)) {
        *fx2 = call_eval_fd(rngi, fdobj2, Lfdobj2);
    }
    
    //apply weight to fx2
    if (call_isempty(wtfd) == false) {
        mxArray *WTD    = call_eval_fd(rngi, wtfd, NULL);
        double *wtd     = mxGetPr(WTD);
        double *fx2d     = mxGetPr(*fx2);
        int M           = mxGetM(WTD);
        int N           = mxGetN(WTD);
        int i, j;
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                fx2d[M*i + j] *= wtd[i];
            }
        }
    }
    
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define TNM_OUT     plhs[0]
    #define S_OUT       plhs[1]
    
    // Macros for input
    #define TNM_IN      prhs[0]
    #define S_IN        prhs[1]
    #define RNGI        prhs[2]
    #define WIDTH       prhs[3]
    #define FDOBJ1      prhs[4]
    #define FDOBJ2      prhs[5]
    #define LFDOBJ1     prhs[6]
    #define LFDOBJ2     prhs[7]
    #define WTFD        prhs[8]
    #define ITER        prhs[9]
    
    double width, *tnm, delta, *x, *rngi, *s, *fx1, *fx2, n;
    int iter, *D, j, k, m, M, index;
    
    TNM_OUT = TNM_IN;
    tnm     = mxGetPr(TNM_OUT);
    tnm[0] *= 2;
    
    width   = mxGetPr(WIDTH)[0];
    delta   = width/tnm[0];
    
    mxArray *X;
    X       = mxCreateDoubleMatrix(1, (int)tnm[0], mxREAL);
    x       = mxGetPr(X);
    
    rngi     = mxGetPr(RNGI);
    
    x[0]    = rngi[0] + delta/2;
    for (j = 1; j < (int)tnm[0]; j++){
        x[j] = x[j-1] + delta;
    }
    
    mxArray *FX1, *FX2;
    
    // evaluate estimated solutions fx1 and fx2 of inner product
    eval_fds(&FX1, &FX2, FDOBJ1, FDOBJ2, LFDOBJ1, LFDOBJ2, WTFD, X);
    
    fx1     = mxGetPr(FX1);
    fx2     = mxGetPr(FX2);
    iter    = (int)mxGetPr(ITER)[0];
    
    // set up s
    S_OUT   = S_IN;
    s       = mxGetPr(S_OUT);
    D       = mxGetDimensions(S_IN);
    M       = mxGetM(FX1);
    
    // assign (iter)th set of solutions
    for (j = 0; j < D[1]; j++) {
        for (k = 0; k < D[2]; k++) {
            n = 0;
            // multiply fx1^T and fx2
            for (m = 0; m < M; m++){
                n += fx1[m + M*j] * fx2[m + M*k];
            }
            index = iter - 1 + D[0]*j + D[0]*D[1]*k;
            s[index] = (s[index - 1] + (n * delta)) / 2;
        }
    }
    
    return;
}