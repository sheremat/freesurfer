function [Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0_m] = mris_gradient_nonlin__read_siemens_coeff(gradfilename)
%MRIS_GRADIENT_NONLIN__READ_SIEMENS_COEFF  spherical harmonics from coeff.grad
%
% [Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0_m] = mris_gradient_nonlin__read_siemens_coeff(gradfilename)
%
%
% Example:
%  
%   mris_gradient_nonlin__setup
%   mris_gradient_nonlin__list_coeff_files
%
%   [Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0_m] = mris_gradient_nonlin__read_siemens_coeff('coeff_AC84.grad');
  
  
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/oct/10
% $Id: mris_gradient_nonlin__read_siemens_coeff.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( isempty(regexp(gradfilename, '\.grad$')) ),
    warning(sprintf('grad file name   [ %s ]   does not end in ".grad" -- may not be Siemens coefficient file...', gradfilename));
  end;


  [fp, message] = fopen(gradfilename, 'rt');

  if ( fp == -1 ),
    error('error opening grad file ''%s'' for reading: "%s"', gradfilename, message);
  end;
  
  disp(sprintf('==> [%s]: reading coefficients from gradient coil file "%s"', mfilename, gradfilename));


  %==--------------------------------------------------------------------==%
  %%% read file header info

  MAX_LINE_LEN = 1024;


  while 1,  % skip blank lines and comments

    commentline = fgets(fp, MAX_LINE_LEN);

    if ( isempty(commentline) ),
      fclose(fp);
      error('could not read the coefficient file %s', gradfilename);
    end;

    if ( strncmp(commentline, '#*] END:', 8) ),
      break;
    end;

  end;

  % skip the next line. (It contains an information about the system type.)
  systemtype_str = fgets(fp, MAX_LINE_LEN);

  disp(sprintf('==> [%s]: reading system type string from coeff.grad file...', mfilename));
  disp(systemtype_str);


  paramline = fgets(fp, MAX_LINE_LEN);

  % check if first paramline contains "win_"; and, if so, parse

  if ( strncmp(paramline, ' win_', 5) ),

    % parse into the four parameters (these don't seem to be used anywhere...)
    [dum, iThreshLo, dum, iThreshUp, dum, iAlgoTyp, dum, iWinDummy] ...
	= strread(paramline, ' %10c%d%13c%d%13c%d%14c%d;');

    % read next line
    paramline = fgets(fp, MAX_LINE_LEN);
  end;

  % only extract radius and ignore rest
  R0_m = sscanf(paramline, '%f', 1);
  warning('returning R0 in units of METERS!');
  
  % read next line, which contains gradient system mode "(0 = typ. tunnel magnet system; 1 = typ. open magnet system)"
  paramline = fgets(fp, MAX_LINE_LEN);
  CoSyMode = sscanf(paramline, '%d', 1);

  % skip the next 5 lines
  for ind = 3:7,
    fgets(fp, MAX_LINE_LEN);
  end;


  %==--------------------------------------------------------------------==%
  %%% begin reading spherical harmonic coefficients

  count = 0;

  while ( ~feof(fp) ),

    coeffline = fgets(fp, MAX_LINE_LEN);

    if ( isempty(deblank(coeffline)) ),
      % occasionally the .grad file ends with several lines of whitespace
      % before the EOF
      break;
    end;
    
    count = count + 1;

    [num(count), A_or_B(count), dum, n(count), dum, m(count), dum, value(count), xyz(count)] ...
        = strread(coeffline, ' %d %1c%1c %d%1c %d%1c     %f     %c');
    %                           1   2  3  4  5  6  7      8      9

  end;

  fclose(fp);


  %==--------------------------------------------------------------------==%
  %%% organize coefficient values

  % initialize coefficient matrices
  [Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z] = deal( zeros(max([(n+1), (m+1)])) );

  % anders's code assumes a fixed number of coeffs (14)
  %[Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z] = deal( zeros(14, 14) );


  coefftype = [A_or_B; xyz].';

  for ind = 1:count,

    coeffstr = coefftype(ind, :);

    switch coeffstr,

     case 'Ax',
      Alpha_x(n(ind)+1, m(ind)+1) = value(ind);

     case 'Ay',
      Alpha_y(n(ind)+1, m(ind)+1) = value(ind);

     case 'Az',
      Alpha_z(n(ind)+1, m(ind)+1) = value(ind);

     case 'Bx',
      Beta_x(n(ind)+1, m(ind)+1) = value(ind);

     case 'By',
      Beta_y(n(ind)+1, m(ind)+1) = value(ind);

     case 'Bz',
      Beta_z(n(ind)+1, m(ind)+1) = value(ind);

     otherwise,
      error('unrecognized coefficient string: "%s"', coeffstr);

    end;  %% switch
  end;  %% for


  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__read_siemens_coeff.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
