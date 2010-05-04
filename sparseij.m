

% [Ir,Jc] = sparseij( S )
% 
%    Returns Ir and Jc components of sparse matrix S.
%
%    Within memory, real sparse matrices are represented by 3 arrays: Pr is
%    an array of non-zero elements in S. Ir is the same length as Pr and
%    contains the row in S occupied by each corresponding element in Pr. 
%    Jc is an array of indices into Pr and Ir, where Jc(j) is the locaion in Pr of
%    the first non-zero element in column j of S. The last element of Jc is the 
%    location of the last element in the last column of S, hence Jc has length equal
%    to the number of columns in S plus 1. Indexing begins from 0 (c-style)
%    rather than 1 (matlab-style).
% 
%    According to this scheme, Ir( Jc(j) + 1 : Jc(j+1) ) + 1 gives the row indices for non-zero
%    elements of column j in S. 
% 
%    See help on sparse matrix representation in the Matlab help documentation for more information. 

% C Kovach 2007