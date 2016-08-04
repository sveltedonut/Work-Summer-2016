//
//  bsplineC.c
//

#include "mex.h"

#define DOFOR(i,to) for(i=0;i<to;i++)
#define jmax 20
#define INDEC(i,j,col) [j+(i)*col]

long bsplvb(double *t, long jhigh, long index, double x,long left, double *biatx)
{
    /* Evaluates b-spline basis functions */
    static long j=0; long jp1,i;
    static double deltal[jmax],deltar[jmax],saved,term;
    /* index=0 start from scratch else j+1 to jout generated*/
    if(!index) {
        j=0;
        biatx[0]=1.;
        if(jhigh<2) return 0;
    }
    for( ; j<jhigh-1;j++) {
        jp1=j+1;
        deltar[j]=t[left+j+1]-x;
        deltal[j]=x-t[left-j];
        saved=0.;
        DOFOR(i,j+1) {
            term=biatx[i]/(deltar[i]+deltal[j-i]);
            biatx[i]=saved+term*deltar[i];
            saved=deltal[j-i]*term;
        }
        biatx[jp1]=saved;
    }
    return(0);
}


/* calculates values and derivatives of a function
 represented by a B-spline uses bsplvb() */

long bsplvd(double *t, long k, double x, long left, double *a,
            double *biatx, long numder) {
    
    long i,il,j,jlow,jpmid,kp1mm,ldummy,m,kp1,mhigh,ideriv;
    double fk,factor,sum;
    
    kp1=k+1; //---------------------------------------------------------------kp = 5
    mhigh = numder < k ? numder:k;
    mhigh = mhigh  > 1 ? mhigh:1;
    /*  call bsplvb  */
    bsplvb(t,kp1-mhigh,0,x,left,a);
    
    mexPrintf("part 1.n.0.0\n");
    
    /*deBoor uses column ordering of FORTRAN, sends biatx not a*/
    if(mhigh==1) {
        DOFOR(j,k) biatx INDEC(0,j,k)=a[j];
        return(0);
    }
    
    ideriv=mhigh;
    for(m=2;m<=mhigh;m++) {
        jpmid=0;
        for(j=ideriv;j<=k;j++) {
            biatx INDEC(j-1,ideriv-1,k)=a[jpmid];
            jpmid++;
        }
        ideriv--;
        bsplvb(t,kp1-ideriv,1,x,left,a);
    }
    DOFOR(j,k) biatx INDEC(j,0,k)=a[j];
    jlow=0;
    DOFOR(i,k) {
        for(j=jlow;j<k;j++)a INDEC(j,i,k)=0.;
        jlow=i;
        a INDEC( i,i,k)=1.;
    }
    
    mexPrintf("part 1.n.0.1\n");
    
    for(m=2;m<=mhigh;m++) {
        kp1mm=kp1-m;
        fk=(double)kp1mm;
        il=left;
        i=k-1;
        DOFOR(ldummy,kp1mm) {
            factor=fk/(t[il+kp1mm]-t[il]);
            DOFOR(j,i+1) a INDEC(i,j,k)=
            factor*(a INDEC(i,j,k)-a INDEC(i-1,j,k));
            il--;i--;
        }
        DOFOR(i,k) {
            sum=0.;
            jlow = i > m-1 ? i:m-1;
            for(j=jlow;j<k;j++)
                sum+=a INDEC(j,i,k) * biatx INDEC(j,m-1,k);
            biatx INDEC(i,m-1,k)=sum;
        }
    }
    return(0);
    
    mexPrintf("part 1.n.0.2\n");
}

long interv(long nbk, double *t, double x, long *left,long *mflag)
{
    /* Determines which interval (between knots) the location x is in  */
    /*   using binary search */
    
    long i,n,top,btm;
    double below,above;
    n=nbk-1; /*index of last value in t[]*/
    if(x<t[0]) {
        *mflag=-1;
        *left=0;
    }
    if(x>=t[n]) {
        *mflag=1;
        *left=n;
    }
    else {
        *mflag=0;
        /* binary search*/
        btm=0;
        top=n;
        while(1) {
            i=(top+btm)>>1;
            below=t[i];
            above=t[i+1];
            if( below<=x) {
                if(x<above)
                {
                    *left=i;
                    return 0;
                }
                /*x>=above*/
                btm=i;
            }
            else top=i;
        }
    }
    return(0);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Macros for output
    #define BSPLINEMAT  plhs[0]
    
    // Macros for input
    #define X           prhs[0]
    #define BREAKS      prhs[1]
    #define NORDER      prhs[2]
    #define NDERIV      prhs[3]
    #define SPARSEWRD   prhs[4]
    
    //void bspline(long n, double *x, long nbreaks, double *breaks, long norder,
             //long nderiv, double *basismat) {
    /*   Matrix of values of Mth derivative of B-splines corresponding to    */
    /*     N args X.                                                         */
    /*                                                                       */
    /*   N          ... number of argument values in X                       */
    /*   X          ... a nondecreasing set of argument values               */
    /*   NBREAKS    ... number of break values (NBASIS=NBREAKS+NORDER-2)     */
    /*   BREAKS     ... a set of strictly increasing break values containing */
    /*                    all X within their range                           */
    /*   NORDER     ... a order of splines (degree + 1)                      */
    /*   NDERIV     ... order of derivative required (< norder)              */
    /*   BSPLINEMAT ... N by NBASIS matrix holding B-spline values           */
    /*                  in this version the matrix is full storage mode.     */
    
    /*  Last modified 1 June 2016  */
    
    // initialize inputs
    long n          = mxGetN(X);
    double *x       = mxGetPr(X);
    long nbreaks    = mxGetN(BREAKS);
    double *breaks  = mxGetPr(BREAKS);
    long norder     = (long)mxGetPr(NORDER)[0];
    long nderiv     = (long)mxGetPr(NDERIV)[0];
    if (nrhs == 5) {
        long sparsewrd = (long)mxGetPr(SPARSEWRD)[0];
    }
    
    long i, j, k, nbasis, left, numder;
    double *knots, *biatx, *temp;
    
    nbasis = nbreaks + norder - 2;
    numder = nderiv + 1;  /*  number of derivatives to compute  */
    
    knots    = (double*)malloc(sizeof(double)*(nbasis+norder));
    biatx    = (double*)malloc(sizeof(double)*norder*numder);
    temp     = (double*)malloc(sizeof(double)*norder);
    
    mexPrintf("part 0.0\n");
    mexPrintf("n = %d, nbreaks = %d, norder = %d, nbasis = %d\n", n, nbreaks, norder, nbasis);
    
    // initialize output
    double *basismat;
    BSPLINEMAT  = mxCreateDoubleMatrix(n, nbasis, mxREAL);
    basismat    = mxGetPr(BSPLINEMAT);
    
    for (i=0; i<n*nbasis; i++) basismat[i] = 0;
    
    for (i=1; i<nbreaks; i++) {
        if (breaks[i]-breaks[i-1] <= 0) {
            //nrerror("The knot sequence breaks should be increasing.");
        }
    }
    
    for (i=1; i<n; i++) {
        if (x[i]-x[i-1] < 0) {
            //nrerror("The point sequence x should be nondecreasing.");
        }
    }

    if (nderiv >= norder) {
        //nrerror("Order of derivative not less than order of spline");
    }

    /*  set up the knot sequence from the break values.  */

    for (k=0;      k<norder;        k++) knots[k] = breaks[0];
    for (k=norder; k<nbasis;        k++) knots[k] = breaks[k-norder+1];
    for (k=nbasis; k<nbasis+norder; k++) knots[k] = breaks[nbreaks-1];
    
    
    mexPrintf("part 0\n");

    /*  loop through the argument values, computing active B-spline values  */

    left = norder - 1;
    for (i=0; i<n; i++) {
        if (x[i] >= breaks[nbreaks-1]) left = nbasis-1;
            else while (knots[left+1] <= x[i]) left++;
        /*  compute the values of active B-splines  at x  */
        bsplvd(knots, norder, x[i], left, temp, biatx, numder);
        
        mexPrintf("part 1.%d.1\n", i);
        
        for (j=0; j<norder; j++) {
            basismat[(left-norder+1+j)*n+i] = biatx[j];
            
            mexPrintf("part 1.%d.2.%d.1\n", i, j);
        }
    }
    
    // free dynamically allocated variables
    free(knots);
    free(biatx);
    free(temp);
}