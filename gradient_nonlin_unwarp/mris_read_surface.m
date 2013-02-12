function [surf_struct, tags, M0_tkr2ras] = mris_read_surface(surf_str)
%MRIS_READ_SURFACE  read FreeSurfer surface and return in *scanner* coords
%
% [surf_struct, tags, M] = mris_read_surface(surf_filename)
%
% wrapper around mris_read_surface_fs.
%
%
% see also MRIS_READ_SURFACE_FS, MRIS_SAVE_SURFACE.
      
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/sep/25
% $Id: mris_read_surface.m,v 1.4 2012/12/22 00:19:59 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  [surf_struct, tags] = mris_read_surface_fs(surf_str);


  % vertex coordinates are stored as "surfaceRAS", so they need to be
  % converted into true RAS (a.k.a., "scannerRAS"). for typical surfaces,
  % this conversion consists only of a translation given by the c_ras
  % position.
%  M0_tkr2ras = mris_read_surface_tkr2ras(tags);
  M0_tkr2ras = mris_show_matrix(tags);
  if ( isempty(M0_tkr2ras) ),
    error('M0_tkr2ras is empty! error in "mris_show_matrix"');
  end;

  surf_struct.surfaceRAS = surf_struct.vertices;
  surf_struct.M0_tkr2ras = M0_tkr2ras;

  % note: this is usually just a shift that removes the centering of the surface at the origin
  surf_struct = mris_transform(surf_struct, M0_tkr2ras);

  surf_struct.tags = tags;
  

  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_read_surface.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
