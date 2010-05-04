


/*
 *  Christopher Kovach, 2007
 *
 *  TrB = blocktrace(SpX,n);
 * 
 *  Computes the sum of n X n blocks in the block-diagonal matrix SpX.
 */


#include "mex.h"
#include "matrix.h"
//#include <string.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
    
    mwSize Nvar = (mwSize) * mxGetPr(prhs[1]);
    mwSize numelInBlock = Nvar*Nvar;


    plhs[0] = mxCreateDoubleMatrix( Nvar, Nvar, (mxComplexity) 0);

    double * xout = mxGetPr(plhs[0]);
    double * xin  = mxGetPr(prhs[0]);
    
        
    for (int i = 0 ; i < mxGetNzmax(prhs[0]) ;  i++)
    {
        
       xout[ i % numelInBlock ] += xin[ i ];
        
    };
    
};

        

