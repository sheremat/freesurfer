function gradfile = mris_gradient_nonlin__pick_coeff_file(gradname)
%MRIS_GRADIENT_NONLIN__PICK_COEFF_FILE
%
% varargout = mris_gradient_nonlin__pick_coeff_file(varargin)
%
%
% see also MRIS_GRADIENT_NONLIN__LIST_COEFF_FILES, MRIS_GRADIENT_NONLIN__READ_SIEMENS_COEFF.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2011/apr/15
% $Id: mris_gradient_nonlin__pick_coeff_file.m,v 1.1 2012/01/18 21:57:33 nicks Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % correspondences between Bays in the Athinoula A. Martinos Center for
  % Biomedical Imaging at MGH to gradient coils (and between magnet systems
  % and gradient coils, e.g., the "timtrio" and the "AS092") is current as
  % of April 2011.

  switch lower(gradname),

   case {'xx00',  'sonata'},
    %
    gradfile = '';

   case {'ac44',  'allegra'},
    %
    gradfile = 'coeff_AC44.grad';

   case {'as05',  'bay2', 'avanto'},
    %
    gradfile = 'coeff_AS05.grad';

   case {'as092', 'bay3', 'timtrio'},
    %
    gradfile = 'coeff_AS092.grad';

   case {'as092', 'bay4', 'timtrio'},
    %
    gradfile = 'coeff_AS092.grad';

   case {'as092', 'bay6', 'timtrio'},
    %
    gradfile = 'coeff_AS092.grad';

   case {'xx00',  'bay7', 'mmr'},
    % 2011/apr --
    gradfile = 'coeff_AC44.grad';

   case {'sc72',  'bay5', '7t'},
    % 2010/nov --
    gradfile = 'coeff_SC72.grad';

   case {'gc99',  'bay8', 'skyra'},
    % 2011/feb -- 2011/aug
    gradfile = 'coeff_GC99.grad';

   case {'ac84'},
    % 2008/jul -- 2010/nov
    gradfile = 'coeff_AC84.grad';

   case {'ac88'},
    % 2007/??? --
    gradfile = 'coeff_AC88.grad';

   case {'as302', 'connectome'},
    % 2011/sep/20 --
    gradfile = 'coeff_AS302.grad';


   otherwise,
    error('unrecognized gradient system name "%s"', gradname);
    mris_gradient_nonlin__list_coeff_files
  end;


  if ( isempty(gradfile) ),
    error('gradient system "%s" recognized but not currently supported', gradname);
  end;
  
  
  
  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/gradient_coil_files/mris_gradient_nonlin__pick_coeff_file.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
