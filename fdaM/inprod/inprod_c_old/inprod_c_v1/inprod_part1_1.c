//
//  inprod_part1_1.c
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
    #define RNGI        plhs[0]
    #define WIDTH       plhs[1]
    #define H_OUT       plhs[2]
    #define S_OUT       plhs[3]
    #define TNM         plhs[4]
    
    // Macros for input
    #define NRNG        prhs[0]
    #define IRNG        prhs[1]
    #define RNGVEC      prhs[2]
    #define JMAX        prhs[3]
    #define NREP1       prhs[4]
    #define NREP2       prhs[5]
    #define FDOBJ1      prhs[6]
    #define FDOBJ2      prhs[7]
    #define LFDOBJ1     prhs[8]
    #define LFDOBJ2     prhs[9]
    #define WTFD        prhs[10]
    
    double *rngi, *width, *rngvec, *h, *s, *tnm, *fx1, *fx2;
    int *iter, jmaxp, irng, nrng, nrep1, nrep2, j, k, M, index;
    
    // retrieve irng and nrng
    irng    = (int)mxGetPr(IRNG)[0];
    nrng    = (int)mxGetPr(NRNG)[0];
    // retrieve rngvec
    rngvec  = mxGetPr(RNGVEC);
    
    // initialize and assign rngi
    RNGI    = mxCreateDoubleMatrix(1, 2, mxREAL);
    rngi    = mxGetPr(RNGI);
    rngi[0] = rngvec[irng-2];
    rngi[1] = rngvec[irng-1];
    
    // adjust rngi
    if (irng > 2) {
        rngi[0] += 1e-10;
    }
    if (irng < nrng){
        rngi[1] -= 1e-10;
    }
    
    // initialize ans calculate width
    WIDTH       = mxCreateDoubleMatrix(1, 1, mxREAL);
    width       = mxGetPr(WIDTH);
    width[0]    = rngi[1] - rngi[0];
    
    // initialize and calculate JMAXP
    jmaxp    = (int)mxGetPr(JMAX)[0] + 1;
    
    
    // retrieve JMAXP, nrep1, nrep2
    nrep1   = (int)mxGetPr(NREP1)[0];
    nrep2   = (int)mxGetPr(NREP2)[0];
    
    // create uninitialized h array
    H_OUT = mxCreateDoubleMatrix(jmaxp, 1, mxREAL);
    //mxSetM(H_OUT, jmaxp_i);
    //mxSetN(H_OUT, 1);
    //mxSetData(H_OUT, mxMalloc(sizeof(double)*jmaxp_i));
    h = mxGetPr(H_OUT);
    
    // assign values to h
    h[0] = 1;
    h[1] = 0.25;
    for (j = 2; j < jmaxp; j++) {
        h[j] = 1;
    }
    
    //initialize solution matrix s
    mwSize D[3]  = {jmaxp, nrep1, nrep2};
    S_OUT           = mxCreateNumericArray(3, D, mxDOUBLE_CLASS, mxREAL);
    s               = mxGetPr(S_OUT);
    
    mxArray *FX1, *FX2;
    
    // evaluate estimated solutions fx1 and fx2 of inner product
    eval_fds(&FX1, &FX2, FDOBJ1, FDOBJ2, LFDOBJ1, LFDOBJ2, WTFD, RNGI);
    
    // initialize tnm
    TNM     = mxCreateDoubleMatrix(1, 1, mxREAL);
    tnm     = mxGetPr(TNM);
    tnm[0]  = 0.5;
    
    // get contents of fx1 and fx2
    fx1     = mxGetPr(FX1);
    fx2     = mxGetPr(FX2);
    
    // assign first set of solutions
    M       = mxGetM(FX1);
    for (j = 0; j < nrep1; j++) {
        for (k = 0; k < nrep2; k++) {
            int m;
            double n = 0;
            // multiply fx1^T and fx2
            for (m = 0; m < M; m++){
                n += fx1[m + M*j] * fx2[m + M*k];
            }
            index = jmaxp*j + jmaxp*nrep1*k;
            s[index] = n*width[0]/2;
            mexPrintf("%f ", s[index]);
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");
    
    return;
}