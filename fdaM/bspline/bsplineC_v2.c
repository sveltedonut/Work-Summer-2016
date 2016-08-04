//
//  bsplineC_v2.c
//  
//
//  Created by Harry Jiang on 2016-06-15.
//
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define BSPLINEMAT  plhs[0]
    
    // Macros for input
    #define X           prhs[0]
    #define BREAKS      prhs[1]
    #define NORDER      prhs[2]
    #define NDERIV      prhs[3]
    
    long i, j, k, r, jj, nxs, nxn, left, leftpr;
    double saved, tr, tl, term, temp;
    
    long nx         = mxGetN(X);
    double *x       = mxGetPr(X);
    long nb         = mxGetN(BREAKS);
    double *breaks  = mxGetPr(BREAKS);
    long norder     = (long)mxGetPr(NORDER)[0];
    long nderiv     = (long)mxGetPr(NDERIV)[0];
    long nd         = nderiv + 1;
    long ns         = nb - 2 + norder;
    
    double *basismat;
    BSPLINEMAT  = mxCreateDoubleMatrix(nx, ns, mxREAL);
    basismat    = mxGetPr(BSPLINEMAT);
    
    for (i = 0; i < nx * ns; i++) basismat[i] = 0;
    
    /*  set up the knot sequence from the break values.  */
    
    double *b       = (double*)malloc(sizeof(double) * nd * norder);
    double *knots   = (double*)malloc(sizeof(double) * (ns + norder));
    
    for (k = 0;      k < norder;        k++) knots[k] = breaks[0];
    for (k = norder; k < ns;            k++) knots[k] = breaks[k-norder+1];
    for (k = ns;     k < ns+norder;     k++) knots[k] = breaks[nb-1];
    
    //for (k = 0; k < ns + norder; k++) mexPrintf("%f ", knots[k]);
    //mexPrintf("\n");
    
    //mexPrintf("nx = %d, nb = %d, nd = %d, norder = %d\n\n", nx, nb, nd, norder);
    
    /*  loop through the argument values, computing active B-spline values  */
    
    left = norder - 1;
    for (i = 0; i < nx; i++) {
        if (x[i] >= breaks[nb - 1]) left = ns - 1;
        else while (knots[left+1] <= x[i]) left++;
        //mexPrintf("left = %d\n", left);
        
        for (j = 0;     j < nd;             j++) b[j] = 1;
        for (j = nd;    j < nd * norder;    j++) b[j] = 0;
        
        nxs = nd - 1;
        
        // bring b up to intended level
        
        for (j = 0; j < norder - nd; j++) {
            saved = 0;
            for (r = 0; r <= j; r++) {
                leftpr              = left + r;
                tr                  = knots[leftpr + 1] - x[i];
                tl                  = x[i] - knots[leftpr - j];
                term                = b[nxs + (r * nd)] / (tr + tl);
                b[nxs + (r * nd)]   = saved + (tr * term);
                saved               = tl * term;
            }
            b[nxs + ((j+1) * nd)] = saved;
        }
        
        /*for (j = 0; j < nd; j++) {
         for (k = 0; k < norder; k++){
         mexPrintf("%f ", b[j + (k * nd)]);
         }
         mexPrintf("\n");
         }
         mexPrintf("\n");*/
        
        // save b-spline values in successive blocks in b
        
        for (jj = 0; jj < nd - 1; jj++) {
            j       = norder - nd + jj;
            saved   = 0;
            nxn     = nxs - 1;
            for (r = 0; r <= j; r++) {
                leftpr              = left + r;
                tr                  = knots[leftpr + 1] - x[i];
                tl                  = x[i] - knots[leftpr - j];
                term                = b[nxs + (r * nd)] / (tr + tl);
                b[nxn + (r * nd)]   = saved + (tr * term);
                saved               = tl * term;
            }
            b[nxn + ((j+1) * nd)] = saved;
            nxs = nxn;
        }
        
        for (jj = nd - 1; jj > 0; jj--) {
            j       = norder - jj;
            nxs = nd - 1;
            for (r = j - 1; r >= 0; r--) {
                leftpr                  = left + r;
                temp                    = (knots[leftpr + 1] - knots[leftpr - j + 1]) / j;
                //mexPrintf("r = %d, l = %d, j = %d\n", leftpr + 1, leftpr - j, j);
                b[nxs + (r * nd)]       = -b[nxs + (r * nd)] / temp;
                b[nxs + ((r + 1) * nd)] = b[nxs + ((r + 1) * nd)] - b[nxs + (r * nd)];
            }
        }
        
        for (j = 0; j < norder; j++) {
            basismat[i + ((left - norder + 1 + j) * nx)] = b[nd - 1 + (j * nd)];
        }
    }
    
    free(b);
    free(knots);
    return;
}