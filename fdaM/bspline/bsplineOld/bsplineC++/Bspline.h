//
//  Bspline.h
//  
//
//  Created by Harry Jiang on 2016-06-13.
//
//

#ifndef Bspline_h
#define Bspline_h

#include <stdio.h>
#include <stdlib.h>

#endif /* Bspline_h */

void bspline(long n, double *x, long nbreaks, double *breaks, long norder,
             long nderiv, double *basismat);
long bsplvb(double *t, long jhigh, long index, double x,long left, double *biatx);
long bsplvd(double *t, long k, double x, long left, double *a,
            double *biatx, long numder);
long interv(long nbk, double *t, double x, long *left,long *mflag);