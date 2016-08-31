//
//  bsplineC.h
//  
//
//  Created by Harry Jiang on 2016-08-30.
//
//

#ifndef bsplineC_h
#define bsplineC_h


#endif /* bsplineC_h */

mxArray* bsplineC(mxArray* X, mxArray* BREAKS, mxArray* NORDER, long nderiv);
mxArray* call_eval_fd(mxArray* RNG, mxArray* BREAK, mxArray* ORDER, mxArray* COEFS, long Lfdobj);
