function [vol, M0_vox2ras, M0_vox2tkr] = mris_read(infile)
%MRIS_READ
%
% varargout = mris_read(varargin)
%
%
% see also MRIS_READ_MGH, MRIS_READ_NII.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/apr/06
% $Id: mris_read.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  [pathstr, namestr, ext] = fileparts(infile);
  
  infiletype = 'unknown';
    
  if ( ~isempty(cell2mat( regexp(infile, {'\.mgz$',    '\.mgh$'}) )) ),
    infiletype = 'mgz';
  end;
  if ( ~isempty(cell2mat( regexp(infile, {'\.nii.gz$', '\.nii$'}) )) ),
    infiletype = 'nii';
  end;
  
  switch infiletype,
   case {'mgz'},
    [vol, M0_vox2ras, M0_vox2tkr] = mris_read_mgh(infile);
   case {'nii'},
    [vol, M0_vox2ras, M0_vox2tkr] = mris_read_nii(infile);
   otherwise,
    error('unsupported input file extension "%s" -- please use "mgh", "mgz", "nii", or "nii.gz"', ext);
  end;
  
  
  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_read.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
