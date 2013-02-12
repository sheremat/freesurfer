function M0_vox2tkr = mris_read_vox2tkr(filename)
%MRIS_READ_VOX2TKR
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/mar/28
% $Id: mris_read_vox2tkr.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  [status, result] = system(sprintf('mri_info --vox2ras-tkr %s', filename));
  M0_vox2tkr = str2num(result);


  % if the orientation of the data is known, matrix is also provided by:
  %    vqm('LIA', [cres rres sres], size(vol), 0);


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_read_vox2tkr.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
