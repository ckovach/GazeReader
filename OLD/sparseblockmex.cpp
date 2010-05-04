

/*
 *  Christopher Kovach, 2007
 *
 *  SpX = sparseblockmex(X,Ir,Jc, nrows);
 * 
 *  Assigns the data in X to sparse matrix SpX with Ir and Jc as given in arguments.
 */


#include "mex.h"
#include "matrix.h"
#include <string.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
    
    if (nrhs != 4) mexErrMsgTxt("SPARSEBLOCKMEX requires 4 arguments");
   
    const mwSize * xsize = mxGetDimensions(prhs[0]);
    const mwSize * JcSize = mxGetDimensions(prhs[2]);
    mwSize ncol;
    if (JcSize[0] > JcSize[1]) ncol = JcSize[0]-1; else ncol = JcSize[1]-1;
    
//    mwSize nblocks = (mwSize) *mxGetPr(prhs[3]);
      mwSize nrows = (mwSize) *mxGetPr(prhs[3]);
    
    nlhs = 1;
    
  
//    plhs[0] = mxCreateSparse( nblocks*xsize[0], ncol,xsize[0]*xsize[1], (mxComplexity) 0);
      plhs[0] = mxCreateSparse( nrows, ncol,xsize[0]*xsize[1], (mxComplexity) 0);
    
    if (plhs[0] == NULL) mexErrMsgTxt("Unable to allocate memory for sparse matrix");
    
    memcpy( mxGetPr(plhs[0]) , mxGetPr(prhs[0]) , xsize[1]*xsize[0]*mxGetElementSize(prhs[0]));

    
    mwIndex * Ir = (mwIndex *)mxMalloc( mxGetN(prhs[1])*mxGetM(prhs[1]) * sizeof(mwIndex) );
    mwIndex * Jc = (mwIndex *)mxMalloc( mxGetN(prhs[2])*mxGetM(prhs[2]) * sizeof(mwIndex) );
        
    double * Irdb = mxGetPr(prhs[1]);
    double * Jcdb = mxGetPr(prhs[2]);
    
    for (int i = 0 ; i< mxGetN(prhs[1])*mxGetM(prhs[1]) ; i++) Ir[i] = (mwIndex) Irdb[i];
    for (int i = 0 ; i< mxGetN(prhs[2])*mxGetM(prhs[2]) ; i++) Jc[i] = (mwIndex) Jcdb[i];
    
    
    mxSetIr(plhs[0], Ir);
    mxSetJc(plhs[0], Jc);
    
}


