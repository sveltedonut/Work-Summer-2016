//
//  inprod_part0_7.c
//  
//
//  Created by Harry Jiang on 2016-05-12.
//
//

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    #define LAST_APPROX 5
    
    // Macros for output
    #define Y_OUT       plhs[0]
    #define CRIT        plhs[1]
    
    // Macros for input
    #define ITER        prhs[0]
    #define S_IN        prhs[1]
    #define H_IN        prhs[2]
    
    double *h, *s, /**xa, *ds, *cs,*/ *y, dy, ho, hp, errval, ssqval, *crit;
    int iter, i, j, k, m, *D, ns, index, w;
    
    iter            = (int)mxGetPr(ITER)[0];
    h               = mxGetPr(H_IN);
    s               = mxGetPr(S_IN);
    D               = mxGetDimensions(S_IN);
    
    double xa[LAST_APPROX];
    double ya[LAST_APPROX];
    double ds[LAST_APPROX];
    double cs[LAST_APPROX];
    
    Y_OUT   = mxCreateDoubleMatrix(D[1], D[2], mxREAL);
    y       = mxGetPr(Y_OUT);
    
    errval  = 0;
    ssqval  = 0;
    CRIT    = mxCreateDoubleMatrix(1, 1, mxREAL);
    crit    = mxGetPr(CRIT);
    
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
    
    if (ssqval > 0) crit[0] = errval/ssqval;
    else            crit[0] = errval;
    
    return;
}