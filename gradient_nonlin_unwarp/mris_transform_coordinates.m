function [B1, B2, B3] = mris_transform_coordinates(A1, A2, A3, Ma2b)
%MRIS_TRANSFORM_COORDINATES  applies affine transformation to list of 3D coordinates
%
% [B1, B2, B3] = mris_transform_coordinates(A1, A2, A3, Ma2b)
%
%
% see also MRIS_TRANSFORM.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/dec/11
% $Id: mris_transform_coordinates.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  %  pA = [A1(:).'; A2(:).'; A3(:).'];
  %  pA = [pA; ones(1, length(pA))];
  %
  %  pB = Ma2b * pA;
  %
  %  B1 = reshape(pB(1, :), size(A1));
  %  B2 = reshape(pB(2, :), size(A2));
  %  B3 = reshape(pB(3, :), size(A3));


  % this implementation (from Anders) is nice since one or more of the input
  % arguments can be scalars (and the others must have the same dimension)

  B1 = Ma2b(1,1)*A1 + Ma2b(1,2)*A2 + Ma2b(1,3)*A3 + Ma2b(1,4);
  B2 = Ma2b(2,1)*A1 + Ma2b(2,2)*A2 + Ma2b(2,3)*A3 + Ma2b(2,4);
  B3 = Ma2b(3,1)*A1 + Ma2b(3,2)*A2 + Ma2b(3,3)*A3 + Ma2b(3,4);


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_transform_coordinates.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
