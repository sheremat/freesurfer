function FS_volgeom = mris_read_surface_volgeom(tags)
%MRIS_READ_SURFACE_VOLGEOM  extracts volume geometry from surface file tags
%
% FS_volgeom = mris_read_surface_volgeom(tags)
%
% for each surface file written to disk, FreeSurfer embeds the direction
% cosines, center position, volume dimensions, and voxel size of the
% reference volume. this can be "mri/filled.mgz" for standard surfaces or
% any other volume (e.g., if the surface was transformed via
% 'mri_surf2surf').
%
%
% see also MRIS_READ_SURFACE.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/jan/16
% $Id: mris_read_surface_volgeom.m,v 1.1 2012/03/11 01:15:33 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % initialize
  FS_volgeom = struct('width', 0, ...
                      'height', 0, ...
                      'depth', 0, ...
                      'xsize', 0, ...
                      'ysize', 0, ...
                      'zsize', 0, ...
                      'x_r', 0, ...
                      'x_a', 0, ...
                      'x_s', 0, ...
                      'y_r', 0, ...
                      'y_a', 0, ...
                      'y_s', 0, ...
                      'z_r', 0, ...
                      'z_a', 0, ...
                      'z_s', 0, ...
                      'c_r', 0, ...
                      'c_a', 0, ...
                      'c_s', 0);


  % extract substrings from tags with the help of regexp
  match = regexp(tags, ...
                 'volume = (?<width>\d+)\s+(?<height>\d+)\s+(?<depth>\d+)', ...
                 'names');
  FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(FS_volgeom, match);

  match = regexp(tags, ...
                 'voxelsize = (?<xsize>\d\.\d+e[\+-]\d\d)\s+(?<ysize>\d\.\d+e[\+-]\d\d)\s+(?<zsize>\d\.\d+e[\+-]\d\d)', ...
                 'names');
  if ( isempty(match) ),
    % old version?
    match = regexp(tags, ...
                   'voxelsize = (?<xsize>\d+\.\d+)\s+(?<ysize>\d+\.\d+)\s+(?<zsize>\d+\.\d+)', ...
                   'names');
  end;
  FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(FS_volgeom, match);

  match = regexp(tags, ...
                 'xras   = (?<x_r>-?\d\.\d+e[\+-]\d\d)\s+(?<x_a>-?\d\.\d+e[\+-]\d\d)\s+(?<x_s>-?\d\.\d+e[\+-]\d\d)', ...
                 'names');
  if ( isempty(match) ),
  match = regexp(tags, ...
                 'xras   = (?<x_r>-?\d+\.\d+)\s+(?<x_a>-?\d+\.\d+)\s+(?<x_s>-?\d+\.\d+)', ...
                 'names');
  end;
  FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(FS_volgeom, match);

  match = regexp(tags, ...
                 'yras   = (?<y_r>-?\d\.\d+e[\+-]\d\d)\s+(?<y_a>-?\d\.\d+e[\+-]\d\d)\s+(?<y_s>-?\d\.\d+e[\+-]\d\d)', ...
                 'names');
  if ( isempty(match) ),
  match = regexp(tags, ...
                 'yras   = (?<y_r>-?\d+\.\d+)\s+(?<y_a>-?\d+\.\d+)\s+(?<y_s>-?\d+\.\d+)', ...
                 'names');
  end;
  FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(FS_volgeom, match);

  match = regexp(tags, ...
                 'zras   = (?<z_r>-?\d\.\d+e[\+-]\d\d)\s+(?<z_a>-?\d\.\d+e[\+-]\d\d)\s+(?<z_s>-?\d\.\d+e[\+-]\d\d)', ...
                 'names');
  if ( isempty(match) ),
  match = regexp(tags, ...
                 'zras   = (?<z_r>-?\d+\.\d+)\s+(?<z_a>-?\d+\.\d+)\s+(?<z_s>-?\d+\.\d+)', ...
                 'names');
  end;
  FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(FS_volgeom, match);

  match = regexp(tags, ...
                 'cras   = (?<c_r>-?\d\.\d+e[\+-]\d\d)\s+(?<c_a>-?\d\.\d+e[\+-]\d\d)\s+(?<c_s>-?\d\.\d+e[\+-]\d\d)', ...
                 'names');
  if ( isempty(match) ),
  match = regexp(tags, ...
                 'cras   = (?<c_r>-?\d+\.\d+)\s+(?<c_a>-?\d+\.\d+)\s+(?<c_s>-?\d+\.\d+)', ...
                 'names');
  end;
  FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(FS_volgeom, match);


  return;



%**************************************************************************%
function FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(FS_volgeom, match)
%MRIS_READ_SURFACE_VOLGEOM__EXTRACT_REGEXP_STRING  helper function
%
% FS_volgeom = mris_read_surface_volgeom__extract_regexp_string(match)
%
%
% see also MRIS_READ_SURFACE_VOLGEOM.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/jan/16
% $Id: mris_read_surface_volgeom.m,v 1.1 2012/03/11 01:15:33 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  names = fieldnames(match);
  for ind = 1:length(names),
    field = names{ind};

    value = str2num(getfield(match, field));

    FS_volgeom = setfield(FS_volgeom, field, value);

  end;


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_read_surface_volgeom.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
