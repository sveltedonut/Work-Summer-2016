//
//  inprodC.c
//  
//
//  Created by Harry Jiang on 2016-05-20.
//
//

#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <math.h>

// call eval.fd from r
SEXP call_eval_fd(SEXP rng, SEXP fdobj, SEXP Lfdobj){
    SEXP statsPackage;
    PROTECT(
        statsPackage = eval(
            lang2( install("getNamespace"), ScalarString(mkChar("stats")) ),
            R_GlobalEnv
    ));
    
    SEXP ans;
    
    if(Lfdobj == NULL){
        PROTECT(
            ans = eval(
                lang3( install("eval.fd"), rng, fdobj),
                statsPackage
        ));
    } else {
        PROTECT(
            ans = eval(
                lang4( install("eval.fd"), rng, fdobj, Lfdobj),
                statsPackage
        ));
    }
    
    UNPROTECT(2);
    return(ans);
}

// evaluate fx1 and fx2
void eval_fds(double**  fx1,      double**  fx2,
              SEXP      fdobj1,   SEXP      fdobj2,
              SEXP      Lfdobj1,  SEXP      Lfdobj2,
              SEXP      wtfd,     SEXP      rng,
              bool      isweighted) {
    
    SEXP FX1, FX2;
    
    // eval_fd for fx1 and fx2 and necessary checks
    FX1 = call_eval_fd(rng, fdobj1, Lfdobj1);
    FX2 = call_eval_fd(rng, fdobj2, Lfdobj2);
    
    double *fx1d     = REAL(FX1);
    double *fx2d     = REAL(FX2);
    
    //apply weight to fx2
    if (isweighted) {
        SEXP WTD        = call_eval_fd(rng, wtfd, NULL);
        SEXP dim        = getAttrib(WTD, R_DimSymbol) ;
        int M           = INTEGER(dim)[0];
        int N           = INTEGER(dim)[1];
        double *wtd     = REAL(WTD);
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

SEXP inprodC(SEXP Rngvec,   SEXP Nrng,
             SEXP Jmax,     SEXP Jmin,
             SEXP Nrep1,    SEXP Nrep2,
             SEXP Fdobj1,   SEXP Fdobj2,
             SEXP Lfdobj1,  SEXP Lfdobj2,
             SEXP Wtfd,     SEXP Isweighted,    SEXP Eps) {
    
    #define LAST_APPROX 5
    
    //initialize variables
    double *rngi, width, *rngvec, *fx1, *fx2, *x, *inprodmat;
    double delta, n, dy, ho, hp, errval, ssqval, crit, eps, w;
    int jmax, jmaxp, irng, nrng, nrep1, nrep2, tnm, iter, i, j, k, m, index, ns, jmin;
    double xa[LAST_APPROX];
    double ya[LAST_APPROX];
    double ds[LAST_APPROX];
    double cs[LAST_APPROX];
    
    // retrieve nrng and rngvec
    rngvec      = REAL(Rngvec);
    nrng        = length(Rngvec);
    
    // initialize rngi
    SEXP Rngi;
    PROTECT(Rngi    = allocVector(REALSXP, 2));
    rngi            = REAL(Rngi);
    
    // get jmax, jmaxp, jmin, nrep1, nrep2, eps
    jmax        = (int)REAL(Jmax)[0];
    jmin        = (int)REAL(Jmin)[0];
    jmaxp       = jmax + 1;
    nrep1       = INTEGER(Nrep1)[0];
    nrep2       = INTEGER(Nrep2)[0];
    eps         = REAL(Eps)[0];
    
    // check weight
    bool isweighted = LOGICAL(Isweighted)[0];
    
    // initialize solution matrices
    double h[jmaxp];
    double s[jmaxp*nrep1*nrep2];
    double y[nrep1*nrep2];
    
    // initialize outputs
    SEXP Inprodmat;
    PROTECT(Inprodmat   = allocMatrix(REALSXP, nrep1, nrep2));
    inprodmat           = REAL(Inprodmat);
    
    for (j = 0; j < nrep1; j++){
        for (k = 0; k < nrep2; k++){
            index            = j + nrep1*k;
            inprodmat[index] = 0;
        }
    }
    
    //loop through each sub-interval
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
        eval_fds(&fx1, &fx2, Fdobj1, Fdobj2, Lfdobj1, Lfdobj2, Wtfd, Rngi, isweighted);
        
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
        // iterate with further subdivisions of integral
        for (iter = 2; iter <= jmax; iter++){
            delta   = width/tnm;
            
            // initialize new evaluation points
            SEXP X;
            PROTECT(X   = allocVector(REALSXP, tnm));
            x           = REAL(X);
            
            // assign new evaluation points
            x[0]    = rngi[0] + delta/2;
            for (j = 1; j < tnm; j++){
                x[j] = x[j-1] + delta;
            }
            
            // evaluate fx1 and fx2
            eval_fds(&fx1, &fx2, Fdobj1, Fdobj2, Lfdobj1, Lfdobj2, Wtfd, X, isweighted);
            
            // assign (iter)th set of approximations
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
                
                // set/reset errval and ssqval
                errval  = 0.0;
                ssqval  = 0.0;
                
                // initialize x-values to evaluate
                for (i = iter - LAST_APPROX; i < iter; i++) {
                    xa[i + LAST_APPROX - iter] = h[i];
                }
                
                // loop through each element of solution matrix
                for (j = 0; j < nrep1; j++) {
                    for (k = 0; k < nrep2; k++) {
                        // hash matrix location to linear location
                        index = j + nrep1*k;
                        
                        // store recent y estimates
                        for (i = iter - LAST_APPROX; i < iter; i++) {
                            ya[i + LAST_APPROX - iter] = s[i + jmaxp*index];
                        }
                        ns = LAST_APPROX - 1;
                        
                        // set difference arrays
                        for (i = 0; i < LAST_APPROX; i++) {
                            ds[i] = ya[i];
                            cs[i] = ya[i];
                        }
                        // set initial values of y
                        y[index] = cs[ns];
                        ns = ns - 1;
                        
                        // polynomial interpolation
                        for (m = 1; m < LAST_APPROX; m++) {
                            //update cs and ds
                            for (i = 0; i < (LAST_APPROX - m); i++) {
                                ho      = xa[i];
                                hp      = xa[i + m];
                                w       = (cs[i+1] - ds[i])/(ho - hp);
                                ds[i]   = hp*w;
                                cs[i]   = ho*w;
                            }
                            
                            // decide whether to fork through highest or lowest estimate
                            if (2*(ns + 1) < (LAST_APPROX - m)) dy = cs[ns + 1];
                            else {
                                dy = ds[ns];
                                ns--;
                            }
                            y[index] += dy;
                        }
                        
                        // check convergence values
                        if (fabs(dy)        > errval) errval = fabs(dy);
                        if (fabs(s[index])  > ssqval) ssqval = fabs(s[jmaxp*index]);
                    }
                }
                
                // test for convergence
                if (ssqval > 0) crit = errval/ssqval;
                else            crit = errval;
                
                if (crit < eps && iter >= jmin) break;
            }
            
            // extrapolation failed; initialize next iteration
            for (j = 0; j < nrep1; j++) {
                for (k = 0; k < nrep2; k++) {
                    index           = jmaxp*j + jmaxp*nrep1*k;
                    s[iter + index] = s[iter + index - 1];
                }
            }
            UNPROTECT(1);
            h[iter] = 0.25 * h[iter - 1];
            tnm    *= 2;
        }
        
        // update approximation to final matrix
        for (j = 0; j < nrep1; j++){
            for (k = 0; k < nrep2; k++){
                index               = j + nrep1*k;
                inprodmat[index]   += y[index];
            }
        }
    }
    
    UNPROTECT(3);
    return (Inprodmat);
}
