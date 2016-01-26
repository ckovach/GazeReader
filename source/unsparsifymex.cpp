

/*
 *  Christopher Kovach, 2007
 *
 * [X,Ir,Jc] = unsparsifymex( SpX ); Returns a full double matrix 
 *  containing the non-zero data in block diagonal sparse matrix SpX 
 *  as well as the Ir and Jc of SpX.
 *
 *
 * ----------- SVN REVISION INFO ------------------
 * $URL$     
 * $Revision$
 * $Date$
 * $Author$
 * ------------------------------------------------
 */


#include "mex.h"
#include "matrix.h"
#include <string.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
    
    if (nrhs != 1) mexErrMsgTxt("UNSPARSIFYMEX requires 1 argument");
    if ( !mxIsSparse(prhs[0]) ) mexErrMsgTxt("First argument must be a sparse matrix.");
   
    mwIndex * Irind = mxGetIr(prhs[0]);
    mwIndex * Jcind = mxGetJc(prhs[0]);
    
    //const mwSize npar = (mwSize)  *mxGetPr(prhs[1]);

    const mwSize npar = Jcind[1]-Jcind[0]; //This of course assumes that npar is constant for all blocks
    
    const mwSize ncol = mxGetN(prhs[0]);

    
    
    plhs[0] =  mxCreateDoubleMatrix(npar,ncol, (mxComplexity) 0);
    
    double * xdata = mxGetPr(prhs[0]);
  
//    mxSetPr( plhs[0],mxGetPr(prhs[0])  );
    memcpy( mxGetPr(plhs[0]) , xdata , ncol*npar*mxGetElementSize(prhs[0]) );
    
   
    if (nlhs > 1) 
    {
        plhs[1] = mxCreateDoubleMatrix( ncol*npar , 1 , (mxComplexity) 0 );
        double * Ir = mxGetPr(plhs[1]);
        for (int i = 0; i < ncol*npar; i++) Ir[i] = (double) Irind[i];
    };

    if (nlhs > 2) 
    {
        plhs[2] = mxCreateDoubleMatrix( ncol+1 , 1 ,  (mxComplexity) 0 );
        double * Jc = mxGetPr(plhs[2]);
        for (int i = 0; i <   1+ncol; i++) Jc[i] = (double) Jcind[i];
    };
    
   
    
}


