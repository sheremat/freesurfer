function varargout = mris_gradient_nonlin__unwarp_volume__batchmode(infile, outfile, gradname, varargin)
%MRIS_GRADIENT_NONLIN__UNWARP_VOLUME__BATCHMODE
%
% mris_gradient_nonlin__unwarp_volume__batchmode(infile, outfile, gradname)
%
% e.g.,
%
%  mris_gradient_nonlin__unwarp_volume__batchmode(...
%      'MPRAGE.mgz', ...
%      'MPRAGE_undis_offline.mgz', ...
%      'coeff_AC84.grad')

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/jan/25
% $Id: mris_gradient_nonlin__unwarp_volume__batchmode.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  ID = '$Id: mris_gradient_nonlin__unwarp_volume__batchmode.m,v 1.2 2012/03/11 01:16:14 jonp Exp $';
  REVISION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  global DEBUG
  if ( ~isempty(DEBUG) ),
    DEBUG = 1
    dbstop if error
  end;


  %==--------------------------------------------------------------------==%

  OPTION__polarity = 'UNDIS';
  OPTION__calc_method = 'direct';
  OPTION__jacobian_correct = 1;
  OPTION__interp = 'cubic';

  OPTION__JacDet_output = 0;
  OPTION__Displacement_output = 0;

  if ( nargin >= 4 && ~isempty(varargin{1}) ),
    OPTION__polarity = varargin{1};
  end;

  if ( nargin >= 5 && ~isempty(varargin{2}) ),
    OPTION__calc_method = varargin{2};
  end;

  if ( nargin >= 6 && ~isempty(varargin{3}) ),
    OPTION__jacobian_correct = str2num(varargin{3});
  end;

  if ( nargin >= 7 && ~isempty(varargin{4}) ),
    OPTION__interp = varargin{4}
  end;

  if ( nargin >= 8 && ~isempty(varargin{5}) ),
    OPTION__JacDet_output = str2num(varargin{5});
  end;

  if ( nargin >= 9 && ~isempty(varargin{6}) ),
    OPTION__Displacement_output = str2num(varargin{6});
  end;


  %==--------------------------------------------------------------------==%

  %pwd
  t0 = timing();

  addpath /usr/local/freesurfer/dev/matlab
  addpath(fileparts(which(mfilename)))
  
  
  %  try,                                           1       2       3       4       5       6
  %                         1       2       3       4       5       6       7       8       9
    disp(sprintf('  %s(''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%d'', ''%s'', ''%d'', ''%d'')', mfilename, infile, outfile, gradname, OPTION__polarity, OPTION__calc_method, OPTION__jacobian_correct, OPTION__interp, OPTION__JacDet_output, OPTION__Displacement_output));

    if ( ~exist(infile, 'file') ),
      error('input volume file "%s" does not exist -- aborting', infile);
    end;

    extind = cell2mat(regexp(infile, {'\.mgz$', '\.mgh$', '\.nii.gz$', '\.nii$'}));
    
    sepind = regexp(infile, '/');
    if ( ~isempty(sepind) ),
      namestr = infile((sepind(end)+1):extind-1);
    else,
      namestr = infile(1:extind-1);
    end;
      
      
    infiletype = 'unknown';

    if ( ~isempty(cell2mat( regexp(infile, {'\.mgz$',    '\.mgh$'}) )) ),
      infiletype = 'mgz';
    end;
    if ( ~isempty(cell2mat( regexp(infile, {'\.nii.gz$', '\.nii$'}) )) ),
      infiletype = 'nii';
    end;

    mr_parms = [];
    switch infiletype,
     case {'mgz'},
      [I, M0rcs2ras, mr_parms] = load_mgh(infile);
     case {'nii'},
      nii = load_nifti(infile);

      I = nii.vol;
      M0rcs2ras = nii.vox2ras;
     otherwise,
      error('unsupported input file extension "%s" -- please use "mgh", "mgz", "nii", or "nii.gz"', infile(extind:end));
    end;


    if ( ndims(I) > 4 ),
      error('input data has %d dimensions -- image data only up to FOUR dimensions is currently supported', ndims(I));
    end;

    if ( isempty(cell2mat( regexp(outfile, {'\.mgz$', '\.mgh$', '\.nii.gz$', '\.nii$'}) )) ),
      error('unsupported output file extension "%s" -- please use "mgh", "mgz", "nii", or "nii.gz"', infile(extind:end));
    end;


    %==--------------------------------------------------------------------==%
    %%% warp the data!

    M1rcs2ras = vox2ras_0to1(M0rcs2ras);

    mris_gradient_nonlin__setup

    if ( ~isempty(regexp(gradname, 'coeff_.*\.grad')) ),
      gradfile = gradname;
    else,
      gradfile = mris_gradient_nonlin__pick_coeff_file(gradname);
    end;

    % JacDet is map defined on original grid, JacDetw is same map
    % interpolated onto warped grid
    % (VR0, VC0, and VS0 are 0-based indices!)
    [Iw, JacDet, JacDetw, VR0, VC0, VS0, VX, VY, VZ] = mris_gradient_nonlin__unwarp_volume(I, ...
                                                  M1rcs2ras, gradfile, ...
                                                  OPTION__calc_method, ...
                                                  OPTION__polarity, ...
                                                  OPTION__jacobian_correct, ...
                                                  OPTION__interp);


    %==--------------------------------------------------------------------==%

    [fp, message] = fopen(outfile, 'w');
    if ( fp == -1 ),
      fclose(fp);
      error(message);
    end;
    fclose(fp);

    pathstr = fileparts(outfile);
    
    extind = cell2mat(regexp(outfile, {'\.mgz$', '\.mgh$', '\.nii.gz$', '\.nii$'}));
    ext = outfile(extind:end);

    % output the auxiliary files with base of input name and path of output name
    outstr = fullfile(pathstr, namestr);
    
    outfiletype = 'unknown';

    if ( ~isempty(cell2mat( regexp(outfile, {'\.mgz$',    '\.mgh$'}) )) ),
      outfiletype = 'mgz';
    end;
    if ( ~isempty(cell2mat( regexp(outfile, {'\.nii.gz$', '\.nii$'}) )) ),
      outfiletype = 'nii';
    end;

    switch outfiletype,
     case {'mgz'},

      save_mgh(Iw, outfile, M0rcs2ras, mr_parms);

      if ( OPTION__JacDet_output ),
        % jacdetinv is inverse jacobian; forward jacobian > 1 indicates a
        % compression of voxel sizes (which is counter-intuitive), whereas
        % inverse jacobian < 1 indicates a compression of voxel sizes---so
        % inverse jacobian is proportional to voxel volume AND inverse
        % warped jacobian when polarity is UNDIS is in true object
        % coordinates.
        save_mgh(JacDet,     regexprep(outfile, ext, ['__jacobian_bias_orig', ext]), M0rcs2ras, mr_parms);
        save_mgh(1./JacDetw, regexprep(outfile, ext, ['__voxel_volumes_warp', ext]), M0rcs2ras, mr_parms);
      end;

      if ( OPTION__Displacement_output ),
%        save_mgh(VR0, regexprep(outfile, ext, ['__VR0', ext]), M0rcs2ras, mr_parms);
%        save_mgh(VC0, regexprep(outfile, ext, ['__VC0', ext]), M0rcs2ras, mr_parms);
%        save_mgh(VS0, regexprep(outfile, ext, ['__VS0', ext]), M0rcs2ras, mr_parms);

        save_mgh(cat(4, VR0, VC0, VS0), [outstr, '__deform_grad_rel', ext], M0rcs2ras, mr_parms);
        save_mgh(cat(4, VX,  VY,  VZ),  [outstr, '__deform_grad_abs', ext], M0rcs2ras, mr_parms);

	% mris_save_m3z(cat(4, VX,  VY,  VZ), ,  [outstr, '__deform_grad_rel.m3z'], M0rcs2ras);

      end;

     case {'nii'},

      if ( strcmp(infiletype, 'nii') ),
        nii.vol = Iw;
        nii.datatype = 16;
        save_nifti(nii, outfile);
      else,
        warning('to write output in NIFTI format input must also be NIFTI; writing MGZ...');
        save_mgh(Iw, fullfile(pathstr, [outstr, '.mgz']), M0rcs2ras);
      end;


      if ( OPTION__JacDet_output ),
        nii.vol = JacDet;
        save_nifti(nii, regexprep(outfile, ext, ['__jacobian_bias_orig', ext]));

        nii.vol = 1./JacDetw;
        save_nifti(nii, regexprep(outfile, ext, ['__voxel_volumes_warp', ext]));
      end;


      if ( OPTION__Displacement_output ),
%        nii.vol = VR0;
%        save_nifti(nii, regexprep(outfile, ext, ['__VR0', ext]));
%
%        nii.vol = VC0;
%        save_nifti(nii, regexprep(outfile, ext, ['__VC0', ext]));
%
%        nii.vol = VS0;
%        save_nifti(nii, regexprep(outfile, ext, ['__VS0', ext]));

        nii4 = nii;
        nii4.dim = [4, size(Iw), 3, 1,1,1];
        nii4.datatype = 16;

        nii4.vol = cat(4, VR0, VC0, VS0);
        save_nifti(nii4, [outstr, '__deform_grad_rel', ext]);

        nii4.vol = cat(4, VX,  VY,  VZ);
        save_nifti(nii4, [outstr, '__deform_grad_abs', ext]);

      end;



     otherwise,
      error('unsupported output file extension "%s" -- please use "mgh", "mgz", "nii", or "nii.gz"', ext);
    end;


  %==--------------------------------------------------------------------==%

%  catch,
%
%    err = lasterror;
%
%    beep
%    fprintf(1, '!!! %s\n\n', err.message);
%    disp(err.stack(1))
%    if ( length(err.stack) > 1 ),
%      disp(err.stack(end))
%    end;
%
%    disp(ID);
%    disp(version);
%
%    if ( DEBUG ),
%      keyboard;
%    end;
%
%  end;


  timing(t0);

  % since called in batch mode, always exit function upon completion
  exit


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__unwarp_volume__batchmode.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
