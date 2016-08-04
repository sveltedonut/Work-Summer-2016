//
//  inprod_part1_3.c
//  
//
//  Created by Harry Jiang on 2016-05-19.
//
//

#include "mex.h"
#include <math.h>

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
    #define LAST_APPROX 5
    
    // Macros for output
    #define Y_OUT       plhs[0]
    
    // Macros for input
    #define S_IN        prhs[0]
    #define H_IN        prhs[1]
    #define RNGI        prhs[2]
    #define WIDTH       prhs[3]
    #define FDOBJ1      prhs[4]
    #define FDOBJ2      prhs[5]
    #define LFDOBJ1     prhs[6]
    #define LFDOBJ2     prhs[7]
    #define WTFD        prhs[8]
    #define EPS         prhs[9]
    #define JMAX        prhs[10]
    #define JMIN        prhs[11]
    
    double *y;
    double *s, *h, *x, *rngi, *fx1, *fx2;
    double width, delta, n, dy, ho, hp, errval, ssqval, crit, eps;
    double xa[LAST_APPROX];
    double ya[LAST_APPROX];
    double ds[LAST_APPROX];
    double cs[LAST_APPROX];
    int jmax, tnm, iter, *D, i, j, k, m, index, ns, w, jmin;
    
    jmax    = (int)mxGetPr(JMAX)[0];
    tnm     = 1;
    width   = mxGetPr(WIDTH)[0];
    rngi    = mxGetPr(RNGI);
    
    eps     = mxGetPr(EPS)[0];
    jmin    = (int)mxGetPr(JMIN)[0];
    
    h       = mxGetPr(H_IN);
    s       = mxGetPr(S_IN);
    D       = mxGetDimensions(S_IN);
    
    Y_OUT   = mxCreateDoubleMatrix(D[1], D[2], mxREAL);
    y       = mxGetPr(Y_OUT);
    
    for (iter = 2; iter <= jmax; iter++){
        delta   = width/tnm;
        
        mxArray *X;
        X       = mxCreateDoubleMatrix(1, tnm, mxREAL);
        x       = mxGetPr(X);
        
        x[0]    = rngi[0] + delta/2;
        for (j = 1; j < tnm; j++){
            x[j] = x[j-1] + delta;
        }
        
        mxArray *FX1, *FX2;
        
        // evaluate estimated solutions fx1 and fx2 of inner product
        eval_fds(&FX1, &FX2, FDOBJ1, FDOBJ2, LFDOBJ1, LFDOBJ2, WTFD, X);
        
        fx1     = mxGetPr(FX1);
        fx2     = mxGetPr(FX2);
        int M   = mxGetM(FX1);
        mexPrintf("M(%d) = %d, tnm = %d\n", iter, M, tnm);
        
        // assign (iter)th set of solutions
        for (j = 0; j < D[1]; j++) {
            for (k = 0; k < D[2]; k++) {
                n = 0;
                // multiply fx1^T and fx2
                for (m = 0; m < tnm; m++){
                    n += fx1[m + tnm*j] * fx2[m + tnm*k];
                }
                index = iter - 1 + D[0]*j + D[0]*D[1]*k;
                s[index] = (s[index - 1] + (n * delta)) / 2;
                mexPrintf("%f ", n);
            }
            mexPrintf("\n");
        }
        mexPrintf("\n");
        
        // Richardson extrapolation
        if (iter >= 5) {
            
            errval  = 0;
            ssqval  = 0;
            
            for (i = iter - LAST_APPROX; i < iter; i++) {
                xa[i + LAST_APPROX - iter] = h[i];
            }
            
            for (j = 0; j < D[1]; j++) {
                for (k = 0; k < D[2]; k++) {
                    // hash matrix location to linear location
                    index = j + D[1]*k;
                    
                    // store recent y estimates
                    for (i = iter - LAST_APPROX; i < iter; i++) {
                        ya[i + 5 - iter] = s[i + D[0]*index];
                    }
                    ns = LAST_APPROX - 1;
                    
                    for (i = 0; i < LAST_APPROX; i++) {
                        ds[i] = ya[i];
                        cs[i] = ya[i];
                    }
                    y[index] = cs[ns];
                    ns--;
                    
                    // Polynomial interpolation
                    for (m = 1; m < LAST_APPROX; m++) {
                        for (i = 0; i < LAST_APPROX - m; i++) {
                            ho      = xa[i];
                            hp      = xa[i + m];
                            w       = (cs[i+1] - ds[i])/(ho - hp);
                            ds[i]   = hp*w;
                            cs[i]   = ho*w;
                        }
                        
                        if (2*ns < LAST_APPROX - m) dy = cs[ns + 1];
                        else {
                            dy = ds[ns];
                            ns--;
                        }
                        y[index] += dy;
                    }
                    
                    // check convergence values
                    if (fabs(dy)        > errval) errval = fabs(dy);
                    if (fabs(s[index])  > ssqval) ssqval = fabs(s[index]);
                }
            }
            
            if (ssqval > 0) crit = errval/ssqval;
            else            crit = errval;
            
            if (crit < eps && iter >= jmin) break;
        }
        
        for (j = 0; j < D[1]; j++) {
            for (k = 0; k < D[2]; k++) {
                index = D[0]*j + D[0]*D[1]*k;
                s[iter + index] = s[iter + index - 1];
            }
        }
        
        h[iter] = 0.25 * h[iter - 1];
        tnm    *= 2;
    }
    
    return;
}