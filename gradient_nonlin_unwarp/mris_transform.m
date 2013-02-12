function surf_struct = mris_transform(surf_struct, M)
%MRIS_TRANSFORM  applies affine transformation to surface vertex coordinates
%
% surf_struct_transformed = mris_transform(surf_struct, M)
%
%
% see also MRIS_TRANSFORM_COORDINATES.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/apr/06
% $Id: mris_transform.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  surf_struct.vertices = M * [surf_struct.vertices, repmat(1, [size(surf_struct.vertices, 1), 1])].';
  surf_struct.vertices =      surf_struct.vertices(1:3,:)';


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_transform.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
