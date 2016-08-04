//
//  inprod_3.c
//  
//
//  Created by Harry Jiang on 2016-05-20.
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
    //mexPrintf("is a fd: %d\n", isa_fjmaxp);
    
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
void eval_fds(double** fx1,             double** fx2,
              const mxArray*  fdobj1,   const mxArray* fdobj2,
              const mxArray*  Lfdobj1,  const mxArray* Lfdobj2,
              const mxArray*  wtfd,     mxArray* rngi) {
    
    mxArray *FX1, *FX2;
    
    // eval_fd for fx1 and fx2 and necessary checks
    if (call_isa_fd(fdobj1)) {
        FX1 = call_eval_fd(rngi, fdobj1, Lfdobj1);
    }
    
    if (call_isa_fd(fdobj2)) {
        FX2 = call_eval_fd(rngi, fdobj2, Lfdobj2);
    }
    double *fx1d     = mxGetPr(FX1);
    double *fx2d     = mxGetPr(FX2);
    
    //apply weight to fx2
    if (call_isempty(wtfd) == false) {
        mxArray *WTD    = call_eval_fd(rngi, wtfd, NULL);
        double *wtd     = mxGetPr(WTD);
        int M           = mxGetM(WTD);
        int N           = mxGetN(WTD);
        int i, j;
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                fx2d[M*i + j] *= wtd[i];
            }
        }
    }
    
    *fx1 = fx1d;
    *fx2 = fx2d;
    
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    #define LAST_APPROX 5
    
    // Macros for output
    #define INPRODMAT   plhs[0]
    
    // Macros for input
    #define NRNG        prhs[0]
    #define RNGVEC      prhs[1]
    #define JMAX        prhs[2]
    #define JMIN        prhs[3]
    #define NREP1       prhs[4]
    #define NREP2       prhs[5]
    #define FDOBJ1      prhs[6]
    #define FDOBJ2      prhs[7]
    #define LFDOBJ1     prhs[8]
    #define LFDOBJ2     prhs[9]
    #define WTFD        prhs[10]
    #define EPS         prhs[11]
    
    double *rngi, width, *rngvec, *fx1, *fx2, *x, *inprodmat;
    double delta, n, dy, ho, hp, errval, ssqval, crit, eps;
    int jmax, jmaxp, irng, nrng, nrep1, nrep2, tnm, iter, *D, i, j, k, m, index, ns, w, jmin;
    double xa[LAST_APPROX];
    double ya[LAST_APPROX];
    double ds[LAST_APPROX];
    double cs[LAST_APPROX];
    
    // retrieve nrngand rngvec
    nrng        = (int)mxGetPr(NRNG)[0];
    rngvec      = mxGetPr(RNGVEC);
    
    // initialize rngi
    mxArray *RNGI;
    RNGI        = mxCreateDoubleMatrix(1, 2, mxREAL);
    rngi        = mxGetPr(RNGI);
    
    // get jmax, jmaxp, jmin, nrep1, nrep2, eps
    jmax        = (int)mxGetPr(JMAX)[0];
    jmin        = (int)mxGetPr(JMIN)[0];
    jmaxp       = jmax + 1;
    nrep1       = (int)mxGetPr(NREP1)[0];
    nrep2       = (int)mxGetPr(NREP2)[0];
    eps         = mxGetPr(EPS)[0];
    
    // initialize solution matrices
    double h[jmaxp];
    double s[jmaxp*nrep1*nrep2];
    double y[nrep1*nrep2];
    
    // initialize outputs
    INPRODMAT   = mxCreateDoubleMatrix(nrep1, nrep2, mxREAL);
    inprodmat   = mxGetPr(INPRODMAT);
    
    for (j = 0; j < nrep1; j++){
        for (k = 0; k < nrep2; k++){
            index            = j + nrep1*k;
            inprodmat[index] = 0;
        }
    }
    
    for (irng = 2; irng <= nrng; irng++){
        // assign rngi
        rngi[0] = rngvec[irng-2];
        rngi[1] = rngvec[irng-1];
        
        // adjust rngi
        if (irng > 2) {
            rngi[0] += 1e-10;
        }
        if (irng < nrng){
            rngi[1] -= 1e-10;
        }
        
        // initialize and calculate width
        width    = rngi[1] - rngi[0];
        
        h[0] = 1;
        h[1] = 0.25;
        
        // evaluate estimated solutions fx1 and fx2 of inner product
        eval_fds(&fx1, &fx2, FDOBJ1, FDOBJ2, LFDOBJ1, LFDOBJ2, WTFD, RNGI);
        
        // assign first set of solutions
        for (j = 0; j < nrep1; j++) {
            for (k = 0; k < nrep2; k++) {
                n = 0.0;
                // multiply fx1^T and fx2
                for (m = 0; m < 2; m++){
                    n += fx1[m + 2*j] * fx2[m + 2*k];
                }
                index    = jmaxp*j + jmaxp*nrep1*k;
                s[index] = n*width/2;
            }
        }
        
        tnm     = 1;
        
        for (iter = 2; iter <= jmax; iter++){
            delta   = width/tnm;
            
            mxArray *X;
            X       = mxCreateDoubleMatrix(1, tnm, mxREAL);
            x       = mxGetPr(X);
            
            x[0]    = rngi[0] + delta/2;
            for (j = 1; j < tnm; j++){
                x[j] = x[j-1] + delta;
            }
            
            // evaluate estimated solutions fx1 and fx2 of inner product
            eval_fds(&fx1, &fx2, FDOBJ1, FDOBJ2, LFDOBJ1, LFDOBJ2, WTFD, X);
            
            // assign (iter)th set of solutions
            for (j = 0; j < nrep1; j++) {
                for (k = 0; k < nrep2; k++) {
                    n = 0.0;
                    // multiply fx1^T and fx2
                    for (m = 0; m < tnm; m++){
                        n += fx1[m + tnm*j] * fx2[m + tnm*k];
                    }
                    index    = iter - 1 + jmaxp*j + jmaxp*nrep1*k;
                    s[index] = (s[index - 1] + (n * delta)) / 2;
                }
            }
            
            // Richardson extrapolation
            if (iter >= 5) {
                
                errval  = 0;
                ssqval  = 0;
                
                for (i = iter - LAST_APPROX; i < iter; i++) {
                    xa[i + LAST_APPROX - iter] = h[i];
                }
                
                for (j = 0; j < nrep1; j++) {
                    for (k = 0; k < nrep2; k++) {
                        // hash matrix location to linear location
                        index = j + nrep1*k;
                        
                        // store recent y estimates
                        for (i = iter - LAST_APPROX; i < iter; i++) {
                            ya[i + 5 - iter] = s[i + jmaxp*index];
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
            
            for (j = 0; j < nrep1; j++) {
                for (k = 0; k < nrep2; k++) {
                    index           = jmaxp*j + jmaxp*nrep1*k;
                    s[iter + index] = s[iter + index - 1];
                }
            }
            
            h[iter] = 0.25 * h[iter - 1];
            tnm    *= 2;
        }
        
        for (j = 0; j < nrep1; j++){
            for (k = 0; k < nrep2; k++){
                index               = j + nrep1*k;
                inprodmat[index]   += y[index];
            }
        }
    }
    
    return;
}
