function varargout = mris_gradient_nonlin__unwarp_surface__batchmode(infile, outfile, gradname, varargin)
%MRIS_GRADIENT_NONLIN__UNWARP_SURFACE__BATCHMODE
%
% mris_gradient_nonlin__unwarp_surface__batchmode(infile, outfile, gradname, regfile, template)
%
% e.g.,
%
%  mris_gradient_nonlin__unwarp_surface__batchmode(...
%      'lh.white', ...
%      'lh.white_DIS', ...
%      'coeff_AC84.grad', ...
%      'register.dat', ...
%      'epi_refframe.nii.gz');

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/apr/01
% $Id: mris_gradient_nonlin__unwarp_surface__batchmode.m,v 1.3 2012/12/22 00:20:16 jonp Exp $
%**************************************************************************%

  ID = '$Id: mris_gradient_nonlin__unwarp_surface__batchmode.m,v 1.3 2012/12/22 00:20:16 jonp Exp $';
  REVISION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG
  if ( isempty(DEBUG) ),
    DEBUG = 0;
  end;

  if ( DEBUG ), dbstop if error; end;

      
  %==--------------------------------------------------------------------==%

  % since for fMRI, typically surfaces are warped *into* distortion, default
  % polarity is "DIS"
  OPTION__polarity = 'DIS';
  OPTION__calc_method = 'direct';
  OPTION__regfile = '';
  OPTION__template = '';
  OPTION__write_area_ratio_overlay = 1;
  OPTION__break = 0;

  if ( nargin >= 4 && ~isempty(varargin{1}) ),
    OPTION__polarity = varargin{1};
  end;

  if ( nargin >= 5 && ~isempty(varargin{2}) ),
    OPTION__calc_method = varargin{2};
  end;

  if ( nargin >= 6 && ~isempty(varargin{3}) ),
    OPTION__regfile = varargin{3};
  end;

  if ( nargin >= 7 && ~isempty(varargin{4}) ),
    OPTION__template = varargin{4};
  end;

  if ( nargin >= 8 && ~isempty(varargin{5}) ),
    OPTION__write_area_ratio_overlay = varargin{5};
  end;

  if ( nargin >= 9 && ~isempty(varargin{6}) ),
    OPTION__break = varargin{6};
  end;


  if ( OPTION__break == 2 ),
    keyboard;
  end;

  pwd
  t0 = timing();

  addpath /usr/local/freesurfer/dev/matlab
  addpath(fileparts(which(mfilename)))

  
  try,

    disp(sprintf('  %s(''%s'', ''%s'', ''%s'', ''%s'', ''%s'')', mfilename, infile, outfile, gradname, OPTION__regfile, OPTION__template));

    if ( ~exist(infile, 'file') ),
      error('input surface file "%s" does not exist -- aborting', infile);
    end;

    [surf_orig, tags, M_surf2RAS] = mris_read_surface(infile);

    if ( ~isempty(OPTION__regfile) ),
      % register.dat is a.k.a. R0_anat2func_tkr, i.e., it is anat2func!
      [regmat, subject, inres, betres] = fmri_readreg(OPTION__regfile);
    else,
      regmat = eye(4);
    end;

    if ( ~isempty(OPTION__template) ),
      [vol, M0_vox2ras, M0_vox2tkr] = mris_read(OPTION__template);
    else,
      M0_vox2tkr = eye(4);
      M0_vox2ras = surf_orig.M0_tkr2ras;
    end;

    %   surface       anat-to-func        func            func
    % (ras-to-tkr) -> (tkr-to-tkr) -> (tkr-to-vox) -> (vox-to-ras)
    surf_orig = mris_transform(surf_orig, M0_vox2ras*inv(M0_vox2tkr)*regmat*inv(surf_orig.M0_tkr2ras));

    
    %==--------------------------------------------------------------------==%
    %%% warp the data!

    mris_gradient_nonlin__setup

    if ( ~isempty(regexp(gradname, 'coeff_.*\.grad')) ),
      gradfile = gradname;
    else,
      gradfile = mris_gradient_nonlin__pick_coeff_file(gradname);
    end;

    if ( OPTION__break == 1 ),
      keyboard;
    end;
    
    % by default, apply distortion (i.e., positive polarity)
    surf_warp = mris_gradient_nonlin__unwarp_surface(surf_orig, gradfile, ...
                                                     OPTION__calc_method, OPTION__polarity);


    %==--------------------------------------------------------------------==%
    %%% calculate area change under deformation

    if ( ~exist('areaPerVertex', 'file') ),
      OPTION__write_area_ratio_overlay = 0;
    end;

    if ( OPTION__write_area_ratio_overlay ),
      surf_orig_area = areaPerVertex(surf_orig);
      surf_warp_area = areaPerVertex(surf_warp);

      surf_area_ratio = surf_warp_area ./ surf_orig_area;
    end;


    %==--------------------------------------------------------------------==%

    [fp, message] = fopen(outfile, 'w');
    if ( fp == -1 ),
      fclose(fp);
      error(message);
    end;
    fclose(fp);

    surf_warp = mris_transform(surf_warp, inv([M0_vox2ras*inv(M0_vox2tkr)*regmat*inv(surf_warp.M0_tkr2ras)]));
    %tags = mris_save_surface__override_tags(tags, OPTION__template, size(vol), voxelsize, M);
    mris_save_surface(outfile, surf_warp, M_surf2RAS, tags);


    %==--------------------------------------------------------------------==%

  catch,

    err = lasterror;

    beep
    fprintf(1, '!!! %s\n\n', err.message);
    disp(err.stack(1))
    if ( length(err.stack) > 1 ),
      disp(err.stack(end))
    end;

    disp(ID);
    disp(version);

    if ( DEBUG ),
      keyboard;
    end;

  end;


  if ( OPTION__write_area_ratio_overlay ),
    save_mgh(surf_area_ratio, [outfile '__area_change_rel.mgz'], eye(4));
  end;


  timing(t0);

  % since called in batch mode, always exit function upon completion
  exit


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__unwarp_surface__batchmode.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
