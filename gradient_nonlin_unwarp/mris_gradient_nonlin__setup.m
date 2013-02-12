function varargout = mris_gradient_nonlin__setup(varargin)
%MRIS_GRADIENT_NONLIN__SETUP  set path to gradient coils coefficients files
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/dec/11
% $Id: mris_gradient_nonlin__setup.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  %if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  mris_toolbox = fileparts(which(mfilename));

  addpath(genpath(mris_toolbox));

  % remove all CVS subdirectories
  path(regexprep(path, ':[^:]*/CVS:', ':'));

  global GRADIENTFILEDIR
  GRADIENTFILEDIR = fullfile(mris_toolbox, 'gradient_coil_files')

  evalin('caller', 'global GRADIENTFILEDIR')


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__setup.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
