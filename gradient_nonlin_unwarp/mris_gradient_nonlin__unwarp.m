function [VXw, VYw, VZw, VJacDet, varargout] = mris_gradient_nonlin__unwarp(VX, VY, VZ, d, gradfile, method, POLARITY)
%MRIS_GRADIENT_NONLIN__UNWARP
%
% [VXw, VYw, VZw, VJacDet] = mris_gradient_nonlin__unwarp(VX, VY, VZ, d, gradfile, method, POLARITY)

% note: if polarity is positive then the coordinates are interpreted as
% ideal or undistorted, and the function returns the locations of these
% points as they would be measured using the specified coil. so,
% counterintuitively, if volume data is being unwarped, then we start with a
% grid of ideal coordinates, find where these locations would fall under the
% imperfect encoding field corresponding to the gradient coil, then stuff
% the intensities found at the measured locations into the ideal locations
% in in the voxel grid (accounting for sub-voxel displacements with
% interpolation).
  
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/dec/13
% $Id: mris_gradient_nonlin__unwarp.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % XYZ (LAI) coordinates of image volume (assume: in *millimeters*)
  V = cat(4, VX, VY, VZ);

  switch method,
   case 'direct',

    gradcoefffile = gradfile;

    % read spherical harmonic coefficients and radius from "coeff.grad" file
    tic;
    [Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0_m] = ...
        mris_gradient_nonlin__read_siemens_coeff(gradcoefffile);
    toc;

    % convert units of radius to millimeters
    R0 = R0_m * 1000;

    % compute the displacements for this set of voxel positions (XYZ coordinates in LAI orientation)
    [DVx, DVy, DVz] = mris_gradient_nonlin__spharm_evaluate(Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0, V);

    % Jacobian determinant is unitless but calculated in terms of
    % displacements and coordinates in LAI orientation (which is right-handed from p.o.v. of **patient**)
    if ( d == 0 ),
      VJacDet = 1;
    else,
      VJacDet = mris_gradient_nonlin__jacobian_determinant(DVx, DVy, DVz, d, VX, VY, VZ);
    end;
    
    %% TODO: grok this; why is any one orientation better than another? siemens uses this one, but why?

    % potential answer: the "magnet coordinates" are defined by the gradient
    % coil and are consistent across scanners. they dictate whether the
    % total field is bigger or smaller than B0 when a positive voltage is
    % applied to the gradient coil. the RAS orientation is specific to the
    % head-first-supine patient positioning, and in the head-first-supine
    % case the magnet coordinates are oriented as LAI. luckily, both RAS and
    % LAI are right-handed, so a rigid rotation is needed to transform
    % between the two, thus the jacobian determinant of this
    % "orientation-preserving" transformation (in the sense of differential
    % geometry) is positive.
    
   case 'lookup',

    gradtablefile = gradfile;

    % displacements are calculated from interpolating values from a
    % previously-generated lookup table for this gradient coil onto the
    % coordinates of the image volume grid
    [DGx, DGy, DGz, GJacDet, G1xyz2rcs, GR, GC, GS] = mris_gradient_nonlin__load_gradtable(gradtablefile);

    % image volume voxels projected into gradient table indices
    [VRg, VCg, VSg] = mris_transform_coordinates(VX, VY, VZ, G1xyz2rcs);

    % values of gradient table evaluated at locations of image voxels 
    % (this interpolation method is fixed at cubic based on internal
    % gradient table file format, so we don't let user pick)
    fprintf(1, '==> [%s]: interpolating...', mfilename);
    tic;
    DVx     = interp3(GR, GC, GS,     DGx, VCg, VRg, VSg, 'cubic');
    DVy     = interp3(GR, GC, GS,     DGy, VCg, VRg, VSg, 'cubic');
    DVz     = interp3(GR, GC, GS,     DGz, VCg, VRg, VSg, 'cubic');

    VJacDet = interp3(GR, GC, GS, GJacDet, VCg, VRg, VSg, 'cubic');
    fprintf(1, 'done! ');
    toc;

   otherwise,
    error('unrecognized displacement calculation method "%s"', method);
  end;  %% switch

  
  % although the inverse mapping is not likely to be simply constructed by
  % subtracting the displacements, if the map is spatially smooth enough
  % this might not be a bad approximation. (the better solution would be to
  % fit a function and invert analytically, or even numerically. it may also
  % be possible to derive it directly from the spherical harmonic
  % coefficients...)
  
  %  - for volumes, a positive POLARITY *removes* the distortion
 
  %  - for surfaces, a positive POLARITY *applies* the distortion
  
  % note that this warping scheme does not lend itself well to removing
  % distortion from surfaces. the removal method starts with the ideal voxel
  % grid then, for each voxel location, it finds the appropriate intensity
  % from the measured image. for surfaces, however, we do not know the ideal
  % vertex locations ahead of time, all we can do is visit each vertex in
  % the measured surface and displace it.
  
  % new locations of image voxels in XYZ (LAI) coordinates
  VXw = VX +  ( POLARITY * DVx );
  VYw = VY +  ( POLARITY * DVy );
  VZw = VZ +  ( POLARITY * DVz );

  % convert jacobian determinant correction if polarity is negative
  if ( sign(POLARITY) == -1 ),
    % cf. VB15 version of MrDistCor2D.cpp, lines 1668 and 1951
    VJacDet = 1./VJacDet;
  end;


  %==--------------------------------------------------------------------==%
  %%% report statistics, and return displacements if requested

  DV = sqrt( DVx(:).^2 + DVy(:).^2 + DVz(:).^2 );

  disp(sprintf('==> [%s]: reporting statistics on displacements...', mfilename));
  disp(sprintf('             max. displacement in 3D = %5.1f mm  (%5.1f mm in x, %5.1f mm in y, %5.1f mm in z )', max(DV),    max(abs(DVx(:))),    max(abs(DVy(:))),    max(abs(DVz(:)))));
  disp(sprintf('             avg. displacement in 3D = %5.1f mm  (%5.1f mm in x, %5.1f mm in y, %5.1f mm in z )', mean(DV),   mean(abs(DVx(:))),   mean(abs(DVy(:))),   mean(abs(DVz(:)))));
  disp(sprintf('             med. displacement in 3D = %5.1f mm  (%5.1f mm in x, %5.1f mm in y, %5.1f mm in z )', median(DV), median(abs(DVx(:))), median(abs(DVy(:))), median(abs(DVz(:)))));


  if ( nargout >= 5 ),
    
    Gx = (DVx + VX)./VX;
    Gy = (DVy + VY)./VY;
    Gz = (DVz + VZ)./VZ;

    varargout{1} = Gx;
    varargout{2} = Gy;
    varargout{3} = Gz;
        
  end;
  

  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__unwarp.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
