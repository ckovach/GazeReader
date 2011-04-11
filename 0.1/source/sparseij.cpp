

/*
 *  Christopher Kovach, 2007
 *
 * [Ir,Jc] = sparseij( S ); Simply returns the Ir and Jc
 *  elements of sparse matrix S. 
 *
 */


#include "mex.h"
#include "matrix.h"
#include <string.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
    
    if (nrhs != 1) mexErrMsgTxt("SPARSEIJ requires 1 argument");
    if ( !mxIsSparse(prhs[0]) ) mexErrMsgTxt("First argument must be a sparse matrix.");
   
    mwIndex * Irind = mxGetIr(prhs[0]);
    mwIndex * Jcind = mxGetJc(prhs[0]);
                          
    const mwSize nnzmax = mxGetNzmax(prhs[0]);
    const mwSize ncol = mxGetN(prhs[0]);
    

    if (nlhs > 0) 
    {
        plhs[0] = mxCreateDoubleMatrix( nnzmax , 1 , (mxComplexity) 0 );
        double * Ir = mxGetPr(plhs[0]);
        for (int i = 0; i < nnzmax; i++) Ir[i] = (double) Irind[i];
    };

    if (nlhs > 1) 
    {
        plhs[1] = mxCreateDoubleMatrix( ncol+1 , 1 ,  (mxComplexity) 0 );
        double * Jc = mxGetPr(plhs[1]);
        for (int i = 0; i <   1+ncol; i++) Jc[i] = (double) Jcind[i];
    };
    
   
    
}


