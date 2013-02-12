function [Iw, varargout] = mris_gradient_nonlin__unwarp_volume(I, M1rcs2ras, gradfile, varargin)
%MRIS_GRADIENT_NONLIN__UNWARP_VOLUME
%
% Iw = mris_gradient_nonlin__unwarp_volume(I, M1rcs2ras, gradfile)
% Iw = mris_gradient_nonlin__unwarp_volume(I, M1rcs2ras, gradfile, METHOD)
% Iw = mris_gradient_nonlin__unwarp_volume(I, M1rcs2ras, gradfile, METHOD, POLARITY)
%
% [Iw, JacDet]  = mris_gradient_nonlin__unwarp_volume(I, M1rcs2ras, gradfile)
% [Iw, JacDetw] = mris_gradient_nonlin__unwarp_volume(I, M1rcs2ras, gradfile)

%
% M1rcs2ras converts 1-based RCS indices into XYZ coordinates in millimeters
% with RAS orientation
%
% polarity is POSITIVE for unwarping, and NEGATIVE for applying a warp; so
%
%  - to unwarp a volume collected with a head-gradient that exhibits strong
%    nonlinearities, the polarity would be POSITIVE.
%
%  - to warp a distortion-free volume (i.e., previously unwarped) into the
%    warped coordinates of a head-gradient, the polarity would be NEGATIVE.
%
% (e.g., "grad_unwarp" and "unwarp_resample.m" use positive polarity)
%
% default polarity is **POSITIVE**: unwarps the volume to remove the distortion

% TODO: account for half-voxel shifts in siemens image coordinates


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/dec/06
% $Id: mris_gradient_nonlin__unwarp_volume.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % by default, calculate displacements from spherical harmonics directly
  % (rather than reading in a previously-calculated lookup table and
  % interpolating to calculate the displacements)
  method = 'direct';
  if ( nargin >= 4 ),
    method = lower(varargin{1});
  end;

  if ( ~[ strcmp(method, 'direct') || strcmp(method, 'lookup') ] ),
    error('unrecognized displacement calculation method "%s": valid options are "direct" or "lookup"', method);
  end;

  if ( ~isempty(regexp(gradfile, '\.gwt$')) ),
    method = 'lookup';
  end;
  if ( ~isempty(regexp(gradfile, '\.grad')) ),
    method = 'direct';
  end;


  % by default. unwarp the volume to remove the distortion
  POLARITY = +1;

  if ( nargin >= 5 ),

    warp_polarity = varargin{2};

    switch upper(warp_polarity),
     case {+1, 'POS', 'FWD', 'UNDIS'},
      % undistort
      POLARITY = +1;
     case {-1, 'NEG', 'INV', 'REV', 'DIS'},
      % "distort"
      POLARITY = -1;
      disp(sprintf('==> [%s]: reverse unwarping polarity selected -- DISTORTING', mfilename));
     otherwise,
      error('invalid warp polarity "%s"; must be ''fwd'' or ''inv''', warp_polarity);
    end;
  end;


  % by default, apply intensity correction based on jacobian determinant of
  % the displacement fields
  FLAG__JACOBIAN_CORRECT = 1;
  if ( nargin >= 6 ),
    FLAG__JACOBIAN_CORRECT = varargin{3};

    if ( FLAG__JACOBIAN_CORRECT ),
      disp(sprintf('==> [%s]: enabling intensity correction based on Jacobian determinant', mfilename));
    else,
      disp(sprintf('==> [%s]: DISABLING intensity correction based on Jacobian determinant', mfilename));
    end;

  end;

  FLAG__CLIP_INTENSITY__SHORT = 0;
  FLAG__CLIP_INTENSITY__FLOAT = 0;

  K_MAX_IMA_PIXEL = 2^12 - 1;  % 4095


  OPTION__interp_method = 'cubic';
  if ( nargin >= 7 ),
    if ( ismember(varargin{4}, {'nearest', 'linear', 'spline', 'cubic'}) ),
      OPTION__interp_method = varargin{4};
    else,
      error('interpolation method must be "nearest", "linear", "spline", or "cubic"');
    end;
  end;
  disp(sprintf('==> [%s]: using %s interpolation', mfilename, upper(OPTION__interp_method)));


  FLAG__gradients_only = 0;
  if ( nargin >= 8 ),
    FLAG__gradients_only = varargin{5};
  end;


  %==--------------------------------------------------------------------==%
  %%%

  % transform RAS-oriented coordinates into LAI-oriented coordinates
  % (magnet coordinates have LAI orientation for the Head-First-Supine position)
  Tras2lai = [-1 0 0 0; 0 +1 0 0; 0 0 -1 0; 0 0 0 1];
  M1rcs2lai = Tras2lai * M1rcs2ras;

  % indices of image volume (one-based)
  [VR, VC, VS] = ndgrid(1:size(I,1), 1:size(I,2), 1:size(I,3));

  % account for half-voxel shift in R and C directions
  % (there is no half-voxel shift in the S direction -- siemens DICOM bug
  % only occurs with in-plane directions!)
  Aedge2centroid = zeros(4);
  Aedge2centroid(1,4) = M1rcs2lai(1,1)/2;
  Aedge2centroid(2,4) = M1rcs2lai(2,2)/2;

  % TODO: check to make sure we are correcting for this shift in the correct
  % direction! e.g., could require a shift to the left or a shift to the
  % right, depending on the in-plane rotation of the slice!
  %M1rcs2lai = M1rcs2lai + Aedge2centroid;


  disp(sprintf('==> [%s]: converting image volume coordinates from RAS to LAI orientation', mfilename));
  [VX, VY, VZ] = mris_transform_coordinates(VR, VC, VS, M1rcs2lai);

  % extract rotation and scaling components of transformation matrix
  % (i.e., ignore translation)
  R1rcs2lai = [   [ M1rcs2lai(1:3,1:3), [0; 0; 0] ];  0 0 0 1   ];

  % since partial derivatives in jacobian matrix are differences that depend
  % on the ordering of the elements of the 3D array, the coordinates may
  % increase in the opposite direction from the array indices, in which case
  % the differential element should be negative. the differentials can be
  % determined by mapping a vector of 1s through the rotation and scaling,
  % where any mirroring will impose a negation
  [dx, dy, dz] = mris_transform_coordinates(1, 1, 1, R1rcs2lai);

  % assemble differentials into one vector to pass to jacobian function
  d = [dx; dy; dz];

  % [dx, dy, dz] = mris_transform_coordinates([1,0,0], [0,1,0], [0,0,1], R0rcs2ras);
  % d = sqrt( dx.^2 + dy.^2 + dz.^2 );


  % (alternatively we could re-write jacobian function to accept M matrix
  % and VR,VC,VS as input instead of VX,VY,VZ)

  if ( FLAG__gradients_only ),

    disp(sprintf('==> [%s]: computing gradients only -- skipping unwarping!', mfilename));
    [null, null, null, null, Gx, Gy, Gz] = ...
        mris_gradient_nonlin__unwarp(VX, VY, VZ, d, gradfile, method, POLARITY);

    Iw = [];

    varargout{1} = Gx;
    varargout{2} = Gy;
    varargout{3} = Gz;

    return;

  else,

    % compute new voxel centroid positions and the jacobian determinant!
    [VXw, VYw, VZw, VJacDetLAI] = mris_gradient_nonlin__unwarp(VX, VY, VZ, d, gradfile, method, POLARITY);

  end;

  % convert back into volume RCS indices
  [VRw, VCw, VSw] = mris_transform_coordinates(VXw, VYw, VZw, inv(M1rcs2lai));

  % calculate displacements (mm units) in original RAS coordinates (**not** LAI!)
  DX = -(VXw - VX);
  DY =  (VYw - VY);
  DZ = -(VZw - VZ);


  % resample the image volume
  % (could also use "interp3(VX, VY, VZ, vol, VXw, VYw, VZw, 'cubic')", but
  % requires more memory)

  % NOTE: interp3 uses MATLAB's "XYZ-convention" for coordinates, which in
  % logical terms corresponds to col-row-slc ordering, thus calls to interp3
  % below are of the form "interp3(I, C, R, S)".

  fprintf(1, '==> [%s]: interpolating...frame', mfilename);
  tic;
  for frame = 1:size(I, 4),
    fprintf(1, ' %d,', frame);
    Iw(:,:,:,frame) = interp3(I(:,:,:,frame), VCw, VRw, VSw, OPTION__interp_method, NaN);
  end;
  fprintf(1, ' done! ');
  toc;


  % assuming that jacobian map is spatially slowly varying compared to the
  % voxel sizes, here we interpolate jacobian calculated on the original
  % grid onto the warped grid
  VJacDetLAIw = interp3(VJacDetLAI, VCw, VRw, VSw, OPTION__interp_method, NaN);


  extrapolated_voxels = find(isnan(Iw));

  disp(sprintf('%s warping resulted in %4.1f%% of voxels shifted outside of image volume', ...
               repmat(' ', length(mfilename) + 7, 1), ...
               length(extrapolated_voxels)/length(I(:)) * 100));


  %% TODO: watch for indices that are outside of volume (voxels are marked
  %% with NaNs by interp3, but can specify another extrapolation value)

  extrapolated_mask = zeros(size(I));
  extrapolated_mask(extrapolated_voxels) = 1;

  % set any NaNs to zero
  Iw(find(isnan(Iw))) = 0;
  VJacDetLAIw(find(isnan(VJacDetLAIw))) = 0;

  if ( FLAG__JACOBIAN_CORRECT ),
    % apply the Jacobian determinant intensity correction
    disp(sprintf('==> [%s]: applying Jacobian determinant correction', mfilename));
    for frame = 1:size(Iw, 4),
      % NOTE: in previous versions this was incorrectly weighted by
      % VJacDetLAI, not VJacDetLAIw!
      Iw(:,:,:,frame) = Iw(:,:,:,frame) .* VJacDetLAIw;
    end;
  end;


  if ( FLAG__CLIP_INTENSITY__SHORT & ~FLAG__CLIP_INTENSITY__FLOAT ),
    % clip all large values to max intensity in a DICOM
    Iw(find( Iw > K_MAX_IMA_PIXEL )) = K_MAX_IMA_PIXEL;
  end;

  if ( ~FLAG__CLIP_INTENSITY__SHORT & FLAG__CLIP_INTENSITY__FLOAT ),

    % heuristic; maybe in future let user change this value, but sooo many inputs already...
    max_intensity = 1.2 * max(vol(:));

    % clip all large values to value slightly higher than max in original volume
    Iw(find( Iw > max_intensity )) = max_intensity;
  end;


  %==--------------------------------------------------------------------==%
  %%% return displacements if requested

  if ( nargout >= 2 ),


    % Jacobian determinant is already expressed in terms of RAS orientation
    varargout{1} = VJacDetLAI;
    varargout{2} = VJacDetLAIw;


    %==------------------------------------------------------------------==%
    %%% prepare displacement maps for FSL tools

    % we need the vox2ras matrix to calculate the voxel sizes
    M0rcs2ras = vox2ras_1to0(M1rcs2ras);

    % remove translation, keep rotation only (so transformed points are
    % distances relative to origin)
    R0rcs2ras = M0rcs2ras; R0rcs2ras(1:3,4) = 0;

    % compute distance traveled along x, y, and z by stepping one voxel in
    % the R, C, and S directions
    [dx, dy, dz] = mris_transform_coordinates([1,0,0], [0,1,0], [0,0,1], R0rcs2ras);
    d = sqrt( dx.^2 + dy.^2 + dz.^2 );

    % NOTE: FSL's "fnirt" tool and associated "applywarp" and "convertwarp"
    % utilities represent displacement or 'warp' fields in an RAS coordinate
    % system with a *radiological* convention for the left-right axis!
    %
    % FSL uses a FOV-based coordinate system in units of *millimeters* (and
    % with the origin at the location of the [0,0,0] voxel) but with the
    % orientation of the coordinate frame established by the vox2ras matrix
    % only, i.e., each value D within the i'th frame of the warp volume file
    % specifies a displacement in the i'th direction by D mm. that is, the
    % direction of the displacements are not absolute but are given by the
    % direction cosines in the vox2ras matrix. the displacement is therefore
    % essentially a voxel or index coordinate system that is scaled by the
    % voxel size.
    %
    % the only caveat is that the tools expect the displacements to be in a
    % *radiological* ordering---even if the image data is not!
    %
    % one can determine whether a file is in neurological ordering with the
    % FSL command "fslorient", e.g.,
    %
    %   fslorient -getorient <infile.nii.gz>
    %
    % if the vox2ras matrix of a image volume (converted from DICOMs with a
    % sane converter like mri_convert, from, say, a product siemens
    % sequence) has a negative determinant then the data is stored in
    % radiological convention whereas if the determinant positive it is
    % neurological. this can be checked on the command line with mri_info:
    %
    %   mri_info -det <infile.nii.gz>
    %
    % so, if the data is neurological, then the sign of the displacement in
    % the ROW dimension of a transverse or coronal image must be negated,
    % e.g., for a neurological ordering if the FOV is 120 mm the matrix is
    % 60 and voxel 13 is to be displaced to the position of voxel 17, the
    % location under the radiological convention is position 103 mm and the
    % relative displacement is -8 mm.

    % so:
    %  - positive determinant => neurological (BAD)
    %  - negative determinant => radiological (GOOD)

    % for readability:
    VR0  = VR-1;
    VC0  = VC-1;
    VS0  = VS-1;

    VRw0 = VRw-1;
    VCw0 = VCw-1;
    VSw0 = VSw-1;

    abs_PR0 = VRw0 * d(1);
    abs_PC0 = VCw0 * d(2);
    abs_PS0 = VSw0 * d(3);

    % NOTE that this flipping is to accommodate FSL, but it assumes that the
    % DICOMs were generated by siemens and follow the standard convention
    % that the L-R axis is the first dimension in transverse and coronal
    % slices. if the data is neurological AND the L-R axis is along a
    % different dimension, then this logic will fail.
    
    % TODO: identify the L-R axis and flip along that. BUT: how will FSL
    % react to this? does it equate the L-R axis with the first dimension,
    % or does it find the axis based on the vox2ras matrix?
    
    % (offline reconstructions from matlab may not obey the strict DICOM
    % convention.)

    % negate relative displacement along ROW dimension if neurological
    rel_PR = d(1) * [VRw-VR];
    rel_PC = d(2) * [VCw-VC];
    rel_PS = d(3) * [VSw-VS];

    % if neurological, flip the voxel indices along the ROW dimension
    if ( det(M0rcs2ras) > 0 ),
      rel_PR = -1 * rel_PR;    % rel_PR = [VR-VRw];
      VR = flipdim(VR, 1);
    end;

    % the absolute displacement is simply the relative displacement plus the
    % voxel locations
    abs_PR = rel_PR + (d(1) * [VR-1]); 
    abs_PC = rel_PC + (d(2) * [VC-1]);
    abs_PS = rel_PS + (d(3) * [VS-1]);

    varargout{3} = rel_PR;
    varargout{4} = rel_PC;
    varargout{5} = rel_PS;

    varargout{6} = abs_PR;
    varargout{7} = abs_PC;
    varargout{8} = abs_PS;


    %==------------------------------------------------------------------==%
    %%% return warped indices (can be passed directly to "interp3.m")

    varargout{09} = VRw;
    varargout{10} = VCw;
    varargout{11} = VSw;


  end;


  return;



  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__unwarp_volume.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
