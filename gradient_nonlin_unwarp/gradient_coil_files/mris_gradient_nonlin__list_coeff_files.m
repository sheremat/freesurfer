function varargout = mris_gradient_nonlin__list_coeff_files(varargin)
%MRIS_GRADIENT_NONLIN__LIST_COEFF_FILES
%

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/dec/11
% $Id: mris_gradient_nonlin__list_coeff_files.m,v 1.1 2012/01/18 21:57:33 nicks Exp $
%**************************************************************************%

  [status, result] = system(sprintf('ls -1 %s/*.grad', fileparts(which(mfilename))));
  if ( status ),
    error(result);
  end;

  fprintf(1, '\nlist of available Siemens gradient coil coefficient files:\n\n');
  disp(result);
  fprintf(1, '\n');
  

  [status, result] = system(sprintf('ls -1 %s/*.gwt', fileparts(which(mfilename))));
  if ( status ),
    error(result);
  end;

  fprintf(1, '\nlist of available Gradient Warp Table (*.gwt) files:\n\n');
  disp(result);
  fprintf(1, '\n');




  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/gradient_coil_files/mris_gradient_nonlin__list_coeff_files.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
