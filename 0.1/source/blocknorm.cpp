


/*
 *  Christopher Kovach, 2007
 *
 *  SpX = blocknorm(E,b);
 * 
 *  Computes the sum of terms within each block and copies to the whole block.
 */


#include "mex.h"
#include "matrix.h"
//#include <string.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
    
    mwSize nblocks;
    const mwSize * bsize = mxGetDimensions(prhs[1]);

    if (bsize[0] > bsize[1]) nblocks = bsize[0]; else nblocks = bsize[1];
    
   // if( mxIsSparse( prhs[0] ) ) mexErrMsgTxt("Input cannot be sparse");
    
    plhs[0] = mxCreateDoubleMatrix( mxGetM(prhs[0]), mxGetN(prhs[0]), (mxComplexity) 0);

    double * xout = mxGetPr(plhs[0]);
    double * xin  = mxGetPr(prhs[0]);
    double * b    = mxGetPr(prhs[1]);
    
    
    double addbuf;
    
    int stindex = 0;
    
    for (int i = 0 ; i < nblocks ;  i++)
    {
        addbuf = 0;
        
        for (int j = 0 ; j < (int) b[i]   ;  j++)
        {
            addbuf += xin[ stindex + j ];
        };
        
        for (int j = 0 ; j < (int) b[i]   ;  j++)
        {
            xout[stindex+j] = addbuf;            
        };
        
        stindex += (int) b[i];
        
    };

    
};

        

