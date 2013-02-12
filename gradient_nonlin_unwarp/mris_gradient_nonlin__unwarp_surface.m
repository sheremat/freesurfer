function [surfStruct_unwarp, varargout] = mris_gradient_nonlin__unwarp_surface(surfStruct, gradfile, varargin)
%MRIS_GRADIENT_NONLIN__UNWARP_SURFACE
%
% surfStruct_unwarp = mris_gradient_nonlin__unwarp_surface(surfStruct, gradfilename, [warp_polarity])
%
% [surfStruct_unwarp, JacDet] = mris_gradient_nonlin__unwarp_surface(surfStruct, gradfilename)
%
%
% polarity is POSITIVE for unwarping, and NEGATIVE for applying a warp; so
%
%  - to unwarp a surface collected with a head-gradient that exhibits strong
%    nonlinearities, the polarity would be POSITIVE.
%
%  - to warp a distortion-free surface (i.e., previously unwarped) into the
%    warped coordinates of a head-gradient, the polarity would be NEGATIVE.
%
% (e.g., "grad_unwarp" and "unwarp_resample.m" use positive polarity)
%
% default polarity is **NEGATIVE**: deforms the surface into a warped field
%
%
% assumes that surface vertices are in scanner coordinates (i.e.,
% translation given by "c_ras" has been applied; see "mris_read_surface.m")
%
% e.g.:
%
%   [surfStruct, M] = mris_read_surface('lh.white');
%   surfStruct_unwarp = mris_gradient_nonlin__unwarp_surface(surfStruct, 'grad_unwarp_tables/ac84.gwv')
%   mris_save_surface(surfStruct_unwarp, M);

% CALLS:
%          mris_gradient_nonlin__load_gradtable

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/oct/10
% $Id: mris_gradient_nonlin__unwarp_surface.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%
  %%% determine polarity of warp: forward or inverse

  % by default, calculate displacements from spherical harmonics directly
  % (rather than reading in a previously-calculated lookup table and
  % interpolating to calculate the displacements)
  method = 'direct';
  if ( nargin >= 3 ),
    method = lower(varargin{1});
  end;

  if ( ~[ strcmp(method, 'direct') || strcmp(method, 'lookup') ] ),
    error('unrecognized displacement calculation method "%s": valid options are "direct" or "lookup"', method);
  end;


  % by default, deform the surface into a warped field
  POLARITY = +1;

  if ( nargin >= 4 ),

    warp_polarity = varargin{2};

    switch upper(warp_polarity),
     case {-1, 'UNDIS', 'REMOVE'},
      % undistort
      POLARITY = -1;
     case {+1, 'DIS', 'APPLY'},
      % "distort"
      POLARITY = +1;
     otherwise,
      error('invalid warp polarity "%s"; must be ''fwd'' or ''inv''', warp_polarity);
    end;
  end;


  %==--------------------------------------------------------------------==%

  
  % to map the surface vertices through a register.dat, remove the
  % transformation stored in the surfStruct "M" matrix, which sets the
  % vertex coordinates in anat_tkr space, apply the foward register.dat,
  % which transforms them to func_tkr space, and apply the inverse
  % M0_func_vox2tkr , and finally the forward M0_func_vox2ras
  
%  M0_func_vox2ras * inv(M0_func_vox2tkr) * reg * inv(surfStruct.M);
  
  %==--------------------------------------------------------------------==%

  % transformation of RAS-oriented coordinates into LAI-oriented coordinates
  % det(Tras2lai) == 1, so maps from "right-handed" to "right-handed"
  Tras2lai = [-1 0 0 0; 0 +1 0 0; 0 0 -1 0; 0 0 0 1];

  % extract coordinates of surface vertices (for readability)
  % (if surface is read using "mris_read_surface", surface vertices are properly translated by c_ras)
  SX = surfStruct.vertices(:,1);
  SY = surfStruct.vertices(:,2);
  SZ = surfStruct.vertices(:,3);

  [VX, VY, VZ] = mris_transform_coordinates(SX, SY, SZ, Tras2lai);

  % NOTE: jacobian defined from surface coordinates does not make sense---needs a voxel size, so best to compute jac for underlying anatomical and paint onto surface separately

%  keyboard
  % compute new coordinates!
  [VXw, VYw, VZw] = mris_gradient_nonlin__unwarp(VX, VY, VZ, 0, gradfile, method, POLARITY);

  % convert back into XYZ with RAS orientation
  [SXw, SYw, SZw] = mris_transform_coordinates(VXw, VYw, VZw, inv(Tras2lai));

  surfStruct_unwarp = surfStruct;
  surfStruct_unwarp.vertices = [ SXw, SYw, SZw ];


  %==--------------------------------------------------------------------==%
  %%% return requested arguments

%  if ( nargout >= 2 ),
%    % return overlay of Jacobian determinant
%    varargout{1} = VJacDetLAI
%  end;

  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__unwarp_surface.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
