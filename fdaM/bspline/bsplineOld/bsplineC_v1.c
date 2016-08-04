//
//  bsplineC.c
//

#include "mex.h"

#define DOFOR(i,to) for(i=0;i<to;i++)
#define jmax 20
#define INDEC(i,j,col) [j+(i)*col]

//original version from code
long bsplvb(double *t, long jhigh, long index, double x,long left, double *biatx)
{
        //Evaluates b-spline basis functions
    static long j; long jp1,i;
    static double deltal[jmax],deltar[jmax],saved,term;
        //index=0 start from scratch else j+1 to jout generated
    if(!index) {
        biatx[0]=1.;
        if(jhigh<2) return 0;
    }
    for(j=0;j<jhigh-1;j++) { //------------------------------------------------ j:               0 -> 1
        jp1=j+1; //----------------------------------------------------------- jp1:             1 -> 2
        deltar[j]=t[left+j+1]-x;
        deltal[j]=x-t[left-j];
        saved=0.;
        DOFOR(i,j+1) { //----------------------------------------------------- i:               0 -> 1
            term=biatx[i]/(deltar[i]+deltal[j-i]);
            biatx[i]=saved+term*deltar[i];
            saved=deltal[j-i]*term;
        }
        biatx[jp1]=saved;
        //mexPrintf("%f ", biatx[jp1]);
    }
    return(0);
}

/*//version from deBoor's book
long bsplvb(double *t, long jhigh, long index, double x,long left, double *biatx)
{
        //Evaluates b-spline basis functions
    static long j; long jp1,i;
    static double deltal[jmax],deltar[jmax],saved,term;
        //index=0 start from scratch else j+1 to jout generated
    if(!index) {
        biatx[0]=1.;
        if(j > jhigh) return 0;
    }
    else{
        for(j=0;j<jhigh-1;j++) { //------------------------------------------------ j:               0 -> 1
            jp1=j+1; //----------------------------------------------------------- jp1:             1 -> 2
            deltar[j]=t[left+j]-x;
            deltal[j]=x-t[left+1-j];
            saved=0;
            DOFOR(i,j) { //----------------------------------------------------- i:               0 -> 1
                term=biatx[i]/(deltar[i]+deltal[jp1-i]);
                biatx[i]=saved+term*deltar[i];
                saved=deltal[jp1-i]*term;
            }
            biatx[jp1]=saved;
            //mexPrintf("%f ", biatx[jp1]);
        }
    }
    return(0);
}*/

/* calculates values and derivatives of a function
 represented by a B-spline uses bsplvb() */

long bsplvd(double *t, long k, double x, long left, double *a, //------------- knots, norder = 4, x[i: 0 -> 10], left, temp, biatx, numder = 2
            double *biatx, long numder) {
    
    // section 0: setup
    long i,il,j,jlow,jpmid,kp1mm,ldummy,m,kp1,mhigh,ideriv, y,z;
    double fk,factor,sum;
    
    kp1=k+1; //--------------------------------------------------------------- kp1 = 3, k = 2
    mhigh = numder < k ? numder:k; //----------------------------------------- mhigh = 2
    mhigh = mhigh  > 1 ? mhigh:1;
    
    //mexPrintf("part 1.n.0.0\n");
    
    // section 1
    /*  call bsplvb  */
    bsplvb(t,kp1-mhigh,0,x,left,biatx);
    //bsplvb(t,k,0,x,left,biatx);
    
    /*DOFOR (i, numder) {
        DOFOR (j, k) {
            mexPrintf("%f ", biatx INDEC(i,j,k));
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");*/
    
    // section 2
    /*deBoor uses column ordering of FORTRAN, sends biatx not a*/
    if(mhigh==1) {
        //DOFOR(j,k) biatx INDEC(0,j,k)=biatx[j];
        return(0);
    }
    /*mexPrintf("biatx:\n");
    DOFOR (y, numder) {
        DOFOR (z, k) {
            mexPrintf("%f ", biatx INDEC(y,z,k));
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");*/

    
    // section 3
    ideriv=mhigh;
    for(m=2;m<=mhigh;m++) {
        jpmid=0;
        for(j=ideriv;j<=k;j++) { //--------------------------------------------- 2 - 2 (runs 1 times

            biatx INDEC(ideriv-1,j-1,k)=biatx[jpmid];
            jpmid++;
            
        }
        ideriv--;
        bsplvb(t,kp1-ideriv,1,x,left,biatx);
    }
    //DOFOR(j,k) biatx INDEC(0,j,k)=a[j]; //------------------------------------ biatx[0 -> 6, 2]     a[0 -> 3]

    
    // section 4
    // set diagonals of a to 1, and top triangle to 0
    jlow=0;
    DOFOR(i,k) {
        for(j=jlow;j<k;j++)a INDEC(i,j,k)=0.;
        jlow=i;
        a INDEC(i,i,k)=1.;
    }
    
    /*mexPrintf("part 1.n.0.4\n");
    
    DOFOR (i, k) {
        DOFOR (j, k) {
            mexPrintf("%f ", a INDEC(j,i,k));
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");*/
    
    // section 5
    for(m=2;m<=mhigh;m++) { //------------------------------------------------ m:           2 -> 2 (runs 1 time)
        kp1mm=kp1-m; //------------------------------------------------------- kp1mm =      3
        fk=(double)kp1mm; //-------------------------------------------------- fk =         3.0
        il=left; //----------------------------------------------------------- il:          3 -> 7 (depends on input)
        i=k-1; //------------------------------------------------------------- i:           3 -> 1
        
        
        
        DOFOR(ldummy,kp1mm) { //---------------------------------------------- ldummy:      0 -> 3 (runs 3 times)
            factor=fk/(t[il+kp1mm]-t[il]); //--------------------------------- factor =     3/(
            DOFOR(j,i) { a INDEC(j,i,k)= //----------------------------------- j:           0 -> 3 (runs 4 times max)           a[15]
                factor*(a INDEC(j,i,k)-a INDEC(j,i-1,k));
                mexPrintf("%d ", j);
                mexPrintf("%d\n", i);
            }
            //mexPrintf("\n");
            il--;i--;
        }
        //mexPrintf("\n");
        
        mexPrintf("a:\n");
        DOFOR (y, k) {
            DOFOR (z, k) {
                mexPrintf("%f ", a INDEC(y,z,k));
            }
            mexPrintf("\n");
        }
        mexPrintf("\n");
        mexPrintf("biatx:\n");
        DOFOR (y, numder) {
            DOFOR (z, k) {
                mexPrintf("%f ", biatx INDEC(y,z,k));
            }
            mexPrintf("\n");
        }
        mexPrintf("\n");
        
        DOFOR(i,k) { //------------------------------------------------------- i:           0 -> 3 (runs 4 times)
            sum=0;
            jlow = i > m-1 ? i :m-1; //---------------------------------------- jlow =       1 -> 3
            for(j=jlow;j<k;j++) {
                sum+= a INDEC(i,j,k) * biatx INDEC(m-1,j,k);
                mexPrintf("%d ", j);
                mexPrintf("%d\n", i);
            }
            biatx INDEC(m-1,i,k)=sum;
            //-------------------------------------- biatx[4 -> 7]
            
            mexPrintf("biatx:\n");
            DOFOR (y, numder) {
                DOFOR (z, k) {
                    mexPrintf("%f ", biatx INDEC(y,z,k));
                }
                mexPrintf("\n");
            }
            mexPrintf("\n");
            //mexPrintf("\n");
        }
        //mexPrintf("\n");
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
    
    long i, j, k, nbasis, left, numder;
    double *knots, *biatx, *temp;
    
    nbasis = nbreaks + norder - 2;
    numder = nderiv + 1;  /*  number of derivatives to compute  */
    
    knots    = (double*)malloc(sizeof(double)*(nbasis+norder)); //double[12]:   [0] -> [11]
    biatx    = (double*)calloc(norder*numder,sizeof(double));   //double[8]:    [0] -> [7]
    temp     = (double*)malloc(sizeof(double)*norder*norder);          //double[16]:    [0] -> [15]
    
    mexPrintf("part 0.0\n");
    mexPrintf("n = %d, nbreaks = %d, norder = %d, nbasis = %d, numder = %d\n", n, nbreaks, norder, nbasis, numder);
    
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
    
    mexPrintf("knots: ");
    for (k = 0;    k<nbasis+norder; k++) mexPrintf("%f ", knots[k]);
    mexPrintf("\n");
    
    
    mexPrintf("part 0\n");

    /*  loop through the argument values, computing active B-spline values  */

    left = norder - 1;
    for (i=0; i<n; i++) {
        if (x[i] >= breaks[nbreaks-1]) left = nbasis-1;
            else while (knots[left+1] <= x[i]) left++; //--------------------- left:            3 -> 7
        mexPrintf("left = %d\n", left);
        /*  compute the values of active B-splines  at x  */
        bsplvd(knots, norder, x[i], left, temp, biatx, numder);
        
        //mexPrintf("part 1.%d.1\n", i);
        
        for (j=0; j<norder; j++) {
            basismat[(left-norder+1+j)*n+i] = biatx[nderiv*norder+j]; //-------------------- basismat[((0 -> 4)+(0 -> 3))*11 + (0 -> 10)]
            mexPrintf("lhs-index = %d, rhs-index = %d", (left-norder+1+j)*n+i, nderiv*norder+j);
            //mexPrintf("part 1.%d.2.%d.1\n", i, j);
        }
        /*DOFOR (j, numder) {
            DOFOR (k, norder) {
                mexPrintf("%f ", biatx[k + norder*j]);
            }
            mexPrintf("\n");
        }
         mexPrintf("\n");*/
        
        DOFOR (j, numder) {
            DOFOR (k, norder) {
                biatx INDEC(j,k,norder) = 0;
            }
        }
    }
    
    // free dynamically allocated variables
    free(knots);
    free(biatx);
    free(temp);
}