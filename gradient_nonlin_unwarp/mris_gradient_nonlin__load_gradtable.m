function [Dx,Dy,Dz, JacDet, G1xyz2rcs, GR, GC, GS, X, Y, Z] = mris_gradient_nonlin__load_gradtable(gradfilename)
%MRIS_GRADIENT_NONLIN__LOAD_GRADTABLE
%
% [Dx,Dy,Dz, JacDet, G1xyz2rcs] = mris_gradient_nonlin__load_gradtable(gradfilename)

% poached from "unwarp_resample.m"

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/jun/24
% $Id: mris_gradient_nonlin__load_gradtable.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( isempty(regexp(gradfilename, '\.gwt$')) ),
    gradfilename = strcat(gradfilename, '.gwt')
  end;

  disp(sprintf('==> [%s]: loading Gradient-nonlinearity Warp Table file "%s"...', mfilename, gradfilename));


  [fp, message] = fopen(gradfilename,'rb','b');

  if ( fp == -1 ),
    disp(sprintf('!!! [%s]: error opening file "%s" for reading', mfilename, gradfilename));
    error(message);
  end;


  % number of grid points along each axis
  ncoords_col = fread(fp, 1, 'int32');
  ncoords_row = fread(fp, 1, 'int32');
  ncoords_slc = fread(fp, 1, 'int32');

  % total number of values
  nv = ncoords_col*ncoords_row*ncoords_slc;

  % extent of grid (in millimeters), given by first and last coordinates [2x1]
  coord_limits_x = fread(fp, 2, 'double');
  coord_limits_y = fread(fp, 2, 'double');
  coord_limits_z = fread(fp, 2, 'double');

  % displacement vectors (in meters)
  dx = fread(fp, nv, 'double').';
  dy = fread(fp, nv, 'double').';
  dz = fread(fp, nv, 'double').';

  % Jacobian determinant vector
  jacDet = fread(fp, nv, 'double')';

  fclose(fp);


  %==--------------------------------------------------------------------==%

  % reshape vectors back into 3D arrays
  Dx = reshape(dx, ncoords_row, ncoords_col, ncoords_slc);
  Dy = reshape(dy, ncoords_row, ncoords_col, ncoords_slc);
  Dz = reshape(dz, ncoords_row, ncoords_col, ncoords_slc);

  JacDet = reshape(jacDet, ncoords_row, ncoords_col, ncoords_slc);


  %==--------------------------------------------------------------------==%

  % length along each axis
  lX = diff(coord_limits_x);
  lY = diff(coord_limits_y);
  lZ = diff(coord_limits_z);
  
  % lX/(ncoords_col-1) is like the size of each grid pixel along X in mm

  % calculate the matrix that maps grid array indices into scanner coordinates
  % (1-based, in units of mm)
  G1rcs2xyz = ...
      [  0,                   lY/(ncoords_row-1),  0,                  coord_limits_y(1)-lY/(ncoords_row-1);
         lX/(ncoords_col-1),  0,                   0,                  coord_limits_x(1)-lX/(ncoords_col-1);
         0,                   0,                   lZ/(ncoords_slc-1), coord_limits_z(1)-lZ/(ncoords_slc-1);
         0, 0, 0, 1];

  % inverse matrix maps scanner coordinates into grid array indices
  % (1-based, in units of mm)
  G1xyz2rcs = inv(G1rcs2xyz);

  % indices into gradient displacement volume
  [GR,GC,GS] = meshgrid(1:ncoords_col, 1:ncoords_row, 1:ncoords_slc);

  
  % reconstitute grid points (in units of mm)
  [X, Y, Z] = meshgrid(...
      linspace(coord_limits_x(1), coord_limits_x(2), ncoords_col), ...
      linspace(coord_limits_y(1), coord_limits_y(2), ncoords_row), ...
      linspace(coord_limits_z(1), coord_limits_z(2), ncoords_slc)  ...
      );


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__load_gradtable.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
