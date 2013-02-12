function [vol, M0_vox2ras, M0_vox2tkr, nii] = mris_read_nii(filename)
%MRIS_READ_NII
%
% [vol, M0_vox2ras, M0_vox2tkr] = mris_read_nii(filename)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/mar/28
% $Id: mris_read_nii.m,v 1.3 2012/12/21 21:55:06 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( ~exist(filename, 'file') ),
    error('file "%s" does not exist -- aborting', filename);
  end;
  
  nii = load_nifti(filename);
  
  vol        = nii.vol; nii.vol = [];
  M0_vox2ras = nii.vox2ras;

  M0_vox2tkr = mris_read_vox2tkr(filename);


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_read_nii.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
