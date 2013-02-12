function mris_save_surface(surf_str, surf_struct, M, varargin)
%MRIS_SAVE_SURFACE
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/sep/25
% $Id: mris_save_surface.m,v 1.3 2012/09/24 05:25:27 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  iM = inv(M);
  surf_struct = mris_transform(surf_struct, iM);

  if ( nargin >= 4 ),
    tags = varargin{1};
  else,
    if ( isfield(surf_struct, 'tags') ),
      tags = getfield(surf_struct, 'tags');
    else
      warning('"tags" field not specified -- using empty string');
      tags = '';
    end;
  end;

  mris_save_surface_fs(surf_struct, surf_str, tags);


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_save_surface.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
