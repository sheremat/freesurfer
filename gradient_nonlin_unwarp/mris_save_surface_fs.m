function varargout = mris_save_surface_fs(surfStruct, surffile, tags)
% MRIS_SAVE_SURFACE_FS  write SGT surface structure to FreeSurfer file
%
% mris_save_surface_fs(surfStruct, surffile, tags)
%
% See also MRIS_SAVE_SURFACE, MRIS_READ_SURFACE, MRIS_READ_SURFACE_FS.

% (based on mukundb's 'readMGHsurf.m', and similar to FreeSurfer's
% "write_surf.m")

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 11/14/2004
% $Id: mris_save_surface_fs.m,v 1.2 2012/12/22 00:19:59 jonp Exp $
%**************************************************************************%

  if ( nargin == 0 ),
    help(mfilename);
    return;
  end;

  [dirname, filename, ext] = fileparts(surffile);

  if ( ~isempty(dirname) & ~exist(dirname, 'dir') ),
    errstr = sprintf('directory "%s" does not exist!', dirname);
    error('==> [%s]: %s\n', mfilename, errstr);
  end;

  if ( exist(surffile, 'file') && getfield(dir(surffile), 'bytes') ),

    backupcmd = sprintf('mv -vf %s %s~', surffile, surffile);
    [status, result] = system(backupcmd);

    warning('SGT:fileClobber', sprintf('file %s exists.', filename));

  end;

  if ( isfield(surfStruct, 'patchfile') ),
    warnstr = sprintf(['input "%s" may be a PATCH struct---writing as' ...
                       ' SURFACE may result in error...'], inputname(1));
    warning('SGT:generic', warnstr);
  end;

  [fid, message] = fopen(surffile, 'wb', 'b');
  if ( fid == -1 ),
    error(message);
  end;


  %==--------------------------------------------------------------------==%

  TRIANGLE_FILE_MAGIC_NUMBER =  16777214;
  QUAD_FILE_MAGIC_NUMBER =  16777215;

  % assume that we are only dealing with triangular meshes for now
  magic = TRIANGLE_FILE_MAGIC_NUMBER;

  fwrite3(fid, magic);

  user_string = getenv('USER');

  date_value = now;
  date_string = sprintf('%s %s %s %s %s', ...
                        datestr(date_value, 'ddd'), ...
                        datestr(date_value, 'mmm'), ...
                        datestr(date_value, 'dd'), ...
                        datestr(date_value, 'HH:MM:SS'), ...
                        datestr(date_value, 'yyyy'));

  fprintf(fid, 'created by %s on %s\n', user_string, date_string);
  fprintf(fid, '\n');

  V = size(surfStruct.vertices,1);
  F = size(surfStruct.faces,   1);

  vertices = surfStruct.vertices;
  faces = surfStruct.faces - 1; % 1-based to 0-based

  vertices = reshape(vertices', [V*3,1]);
  faces = reshape(faces', [F*3,1]);

  fwrite(fid, V, 'int32');
  fwrite(fid, F, 'int32');

  fwrite(fid, vertices, 'float32');
  fwrite(fid, faces, 'int32');


  % TODO: parse the tags and retain only the volume geometry info and, maybe,
  % append info to command line history. example of tag information stored in
  % surface files can be found in the routine "MRISreadOverAlloc()" defined in
  % 'utils/mrisurf.c'. also see 'utils/tags.c' and 'include/tags.h' for more
  % information.
  fwrite(fid, tags, 'uchar');


  status = fclose(fid);

  fprintf(1, '\n');
  dstr = sprintf('file "%s%s" written successfully!', filename, ext);
  disp(sprintf('==> [%s]: %s\n', mfilename, dstr));

  if ( nargout > 1 ),
    varargout{1} = status;
  end;

  return;


%**************************************************************************%
function fwrite3(fid, val)
% write a 3 byte integer out of a file
% (local copy of ${FREESURFER_HOME}/matlab/fwrite3.m)

  b1 = bitand(bitshift(val, -16), 255) ;
  b2 = bitand(bitshift(val, -8), 255) ;
  b3 = bitand(val, 255) ;
  fwrite(fid, b1, 'uchar') ;
  fwrite(fid, b2, 'uchar') ;
  fwrite(fid, b3, 'uchar') ;

  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_save_surface_fs.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: