//
//  bsplineC.c
//  
//
//  Created by Harry Jiang on 2016-06-23.
//
//

#include <R.h>
#include <Rinternals.h>

SEXP bsplineC(SEXP X, SEXP Breaks, SEXP Norder, SEXP Nderiv) {
    // Initialize variables
    long i, j, k, r, jj, nd, ns, nxs, nxn, left, leftpr;
    double saved, tr, tl, term, temp;
    
    // Initialize inputs and constants
    long nx         = length(X);                    // number of evaluation points
    double *x       = REAL(X);                      // array of evaluation points
    long nb         = length(Breaks);               // number of breakpoints
    double *breaks  = REAL(Breaks);                 // array of breakpoints
    long norder     = (long)REAL(Norder)[0];     // order of b-spline
    long nderiv     = (long)REAL(Nderiv)[0];     // derivative requested
    
    nd         = nderiv + 1;
    ns         = nb - 2 + norder;                   // number of splines to evaluate
    
    // Initialize output
    double *basismat;
    SEXP Basismat;
    PROTECT(Basismat    = allocMatrix(REALSXP, nx, ns)); // initialize nx by ns matrix
    basismat            = REAL(Basismat);
    for (i = 0; i < nx * ns; i++) basismat[i] = 0;
    
    
    // Set up the knot sequence from the break values.
    double *b       = (double*)malloc(sizeof(double) * nd * norder);    //Store b-spline values
    double *knots   = (double*)malloc(sizeof(double) * (ns + norder));  //Knots list
    
    for (k = 0;      k < norder;        k++) knots[k] = breaks[0];
    for (k = norder; k < ns;            k++) knots[k] = breaks[k-norder+1];
    for (k = ns;     k < ns+norder;     k++) knots[k] = breaks[nb-1];
    
    // Loop through the argument values, computing active B-spline values
    
    left = norder - 1;
    for (i = 0; i < nx; i++) {
        
        // Find value of left
        if (x[i] >= breaks[nb - 1]) left = ns - 1;
        else while (knots[left+1] <= x[i]) left++;
        
        for (j = 0;     j < nd;             j++) b[j] = 1;
        for (j = nd;    j < nd * norder;    j++) b[j] = 0;
        
        nxs = nd - 1;
        
        // Bring b up to intended level
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
        
        // Save b-spline values b
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
        
        // Now use the fact that derivative values can be obtained by differencing:
        for (jj = nd - 1; jj > 0; jj--) {
            j       = norder - jj;
            nxs = nd - 1;
            for (r = j - 1; r >= 0; r--) {
                leftpr                  = left + r;
                temp                    = (knots[leftpr + 1] - knots[leftpr - j + 1]) / j;
                b[nxs + (r * nd)]       = -b[nxs + (r * nd)] / temp;
                b[nxs + ((r + 1) * nd)] = b[nxs + ((r + 1) * nd)] - b[nxs + (r * nd)];
            }
        }
        
        // Assign last row of b to ith row of basismat
        for (j = 0; j < norder; j++) {
            basismat[i + ((left - norder + 1 + j) * nx)] = b[nd - 1 + (j * nd)];
        }
    }
    
    // Free dynamically allocated variables
    free(b);
    free(knots);
    UNPROTECT(1);
    return(Basismat);
}