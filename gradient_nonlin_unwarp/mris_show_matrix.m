function M = mris_show_matrix(tags)
%MRIS_SHOW_MATRIX  read surface coordinate matrix from a FreeSurfer surface
%
%
% /space/padkeemao/1/users/jonp/lwlab/PROJECTS/VISUOTOPY/mris_toolbox/mris_show_matrix
%
% /space/padkeemao/1/users/jonp/lwlab/PROJECTS/VISUOTOPY/mris_toolbox/show_matrix/mris_show_matrix.c


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/apr/16
% $Id: mris_show_matrix.m,v 1.1 2012/03/11 01:15:33 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%


  FS_volgeom = mris_read_surface_volgeom(tags);

  m = MRIxfmCRS2XYZ(FS_volgeom, 0);

  tmp = FS_volgeom;

  tmp.x_r = -1;
  tmp.y_r =  0;
  tmp.z_r =  0;
  tmp.c_r = 0.0;
  tmp.x_a =  0;
  tmp.y_a =  0;
  tmp.z_a =  1;
  tmp.c_a = 0.0;
  tmp.x_s =  0;
  tmp.y_s = -1;
  tmp.z_s =  0;
  tmp.c_s = 0.0;

  K = MRIxfmCRS2XYZ(tmp,0);

  M = m * inv(K);


  return;



%**************************************************************************%
function m = MRIxfmCRS2XYZ(mri, base)
%MRIxfmCRS2XYZ  poached from FreeSurfer's "mri.c"
%
% M = MRIxfmCRS2XYZ(FS_volgeom)
%
%
% see also MRIS_SHOW_MATRIX.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/jan/16
% $Id: mris_show_matrix.m,v 1.1 2012/03/11 01:15:33 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  m = zeros(4);

  % direction cosine between columns scaled by distance between colums
  m(1, 1) = mri.x_r * mri.xsize;
  m(2, 1) = mri.x_a * mri.xsize;
  m(3, 1) = mri.x_s * mri.xsize;

  % direction cosine between rows scaled by distance between rows
  m(1, 2) = mri.y_r * mri.ysize;
  m(2, 2) = mri.y_a * mri.ysize;
  m(3, 2) = mri.y_s * mri.ysize;

  % direction cosine between slices scaled by distance between slices
  m(1, 3) = mri.z_r * mri.zsize;
  m(2, 3) = mri.z_a * mri.zsize;
  m(3, 3) = mri.z_s * mri.zsize;

  % Preset the offsets to 0
  m(1, 4) = 0.0;
  m(2, 4) = 0.0;
  m(3, 4) = 0.0;

  % Last row of matrix
  m(4, 1) = 0.0;
  m(4, 2) = 0.0;
  m(4, 3) = 0.0;
  m(4, 4) = 1.0;

  % At this point, m = Mdc * D

  % Col, Row, Slice at the Center of the Volume
  Pcrs = zeros(4);
  Pcrs(1, 1) = mri.width/2.0  + base;
  Pcrs(2, 1) = mri.height/2.0 + base;
  Pcrs(3, 1) = mri.depth/2.0  + base;
  Pcrs(4, 1) = 1.0;

  % XYZ offset the first Col, Row, and Slice from Center
  % PxyzOffset = Mdc*D*PcrsCenter
  PxyzOffset = m * Pcrs;

  % XYZ at the Center of the Volume is mri.c_r, c_a, c_s

  % The location of the center of the voxel at CRS = (0,0,0)
  m(1, 4) = mri.c_r - PxyzOffset(1, 1);
  m(2, 4) = mri.c_a - PxyzOffset(2, 1);
  m(3, 4) = mri.c_s - PxyzOffset(3, 1);


  return;

  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_show_matrix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
