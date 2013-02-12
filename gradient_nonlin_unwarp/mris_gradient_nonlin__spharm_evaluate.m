function [Dx, Dy, Dz, X, Y, Z] = mris_gradient_nonlin__spharm_evaluate(Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0, user_coords, resolution)
%MRIS_GRADIENT_NONLIN__DISPLACEMENTS
%
%
%  [Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0] = mris_gradient_nonlin__read_siemens_coeff('coeff.grad');
%  [Dx, Dy, Dz, X, Y, Z] = mris_gradient_nonlin__spharm_evaluate(Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0, 1.2*R0, 'med');
%  JacDet = mris_gradient_nonlin__jacobian_determinant(Dx, Dy, Dz, X, Y, Z);
%
%  mris_gradient_nonlin__save_gradtable('coeff', Dx, Dy, Dz, JacDet, X, Y, Z);


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/oct/10
% $Id: mris_gradient_nonlin__spharm_evaluate.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % note: all units are in millimeters


  if ( exist('resolution', 'var') ),
    if ( isstr(resolution) ),
      switch lower(resolution(1:2)),
       case 'lo',         % low resolution
        resolution_mm = 10;

       case {'me', 'in'}, % medium or intermediate resolution
        resolution_mm =  5;

       case 'hi',         % high resolution
        resolution_mm =  1;

       otherwise,
        error('resolution string "%s" unknown -- options: "low", "med", or "high"', resolution);
      end;
    else,
      resolution_mm = resolution;
    end;
  end;



  if (     isscalar( user_coords ) ),

    disp(sprintf('==> [%s]: generating coordinates from -%2.2f mm to +%2.2f mm', mfilename, user_coords, user_coords));
    coords = -user_coords : resolution_mm : +user_coords;
    [X, Y, Z] = ndgrid(coords, coords, coords);

  elseif ( isvector( user_coords ) && length(user_coords) == 3 ),

    disp(sprintf('==> [%s]: generating X coordinates from -%2.2f mm to +%2.2f mm', mfilename, user_coords(1), user_coords(1)));
    disp(sprintf('==> [%s]: generating Y coordinates from -%2.2f mm to +%2.2f mm', mfilename, user_coords(2), user_coords(2)));
    disp(sprintf('==> [%s]: generating Z coordinates from -%2.2f mm to +%2.2f mm', mfilename, user_coords(3), user_coords(3)));

    coords_X = -user_coords(1) : resolution_mm : +user_coords(1);
    coords_Y = -user_coords(2) : resolution_mm : +user_coords(2);
    coords_Z = -user_coords(3) : resolution_mm : +user_coords(3);

    [X, Y, Z] = ndgrid(coords_X, coords_Y, coords_Z);

  elseif ( ndims(user_coords) == 2 ),  % it is a matrix!

    if ( size(user_coords, 2) == 3 ),
      user_coords = user_coords.';
    end;

    coords_X = user_coords(:, 1);
    coords_Y = user_coords(:, 2);
    coords_Z = user_coords(:, 3);

    [X, Y, Z] = ndgrid(coords_X, coords_Y, coords_Z);

  elseif ( ndims(user_coords) == 4 ),  % it is a stack of 3D arrays!

    disp(sprintf('==> [%s]: using user-supplied coordinates', mfilename));

    X = user_coords(:,:,:,1);
    Y = user_coords(:,:,:,2);
    Z = user_coords(:,:,:,3);

  else,

    error('specified coordinates are in an unknown format');

  end;

  x = X(:).';
  y = Y(:).';
  z = Z(:).';

  disp(sprintf('==> [%s]: calculating displacements (in mm) using spherical harmonic coefficients...', mfilename));
  t0 = timing;
  disp('  computing displacements along x');
  bx = mris_gradient_nonlin__siemens_B(Alpha_x, Beta_x, x, y, z, R0);
  disp('  computing displacements along y');
  by = mris_gradient_nonlin__siemens_B(Alpha_y, Beta_y, x, y, z, R0);
  disp('  computing displacements along z');
  bz = mris_gradient_nonlin__siemens_B(Alpha_z, Beta_z, x, y, z, R0);
  timing(t0);

  Bx = reshape(bx, size(X));
  By = reshape(by, size(Y));
  Bz = reshape(bz, size(Z));

  Dx = Bx * R0;
  Dy = By * R0;
  Dz = Bz * R0;


  return;


%**************************************************************************%
function B = mris_gradient_nonlin__siemens_B(Alpha, Beta, X, Y, Z, R0, Gref)
%MRIS_GRADIENT_NONLIN__SIEMENS_B  calculate displacement field from Siemens's coefficients

  nmax = size(Alpha,1)-1;

  % hack to avoid singularities at origin (R==0)
  X = X+0.0001;

  % convert to spherical coordinates
  R = sqrt(X.^2+Y.^2+Z.^2);
  Theta = acos(Z./R);
  Phi = atan2(Y./R,X./R);

  if ( 1 ),

    % evaluate the Legendre polynomial (using Siemens's normalization)
    B = 0;
    for n = 0:nmax,
      P = mris_gradient_nonlin__siemens_legendre(n,cos(Theta));
      F = (R/R0).^n;
      for m = 0:n,
        F2 = Alpha(n+1,m+1)*cos(m*Phi)+Beta(n+1,m+1)*sin(m*Phi);
        B = B+F.*P(m+1,:).*F2;
      end;
    end;

  else,
%
%       % set up vandermonde matrix and coeffs vector
%       vdm = [];
%       coeffs = [];
%       counter = 0;
%       for n = 0:nmax
%         P = siemens_legendre(n,cos(Theta));
%         F = (R/R0).^n;
%         for m = 0:n
%
%           if( abs(Alpha(n+1,m+1)) > sqrt(eps) )
%             counter = counter + 1;
%             vdm(:, counter) = F.*cos(m*Phi).*P(m+1,:)';
%             coeffs(counter) = Alpha(n+1,m+1);
%           end
%
%           if( abs(Beta(n+1,m+1)) > sqrt(eps) )
%             counter = counter + 1;
%             vdm(:, counter) = F.*sin(m*Phi).*P(m+1,:)';
%             coeffs(counter) = Beta(n+1,m+1);
%           end
%
%         end
%       end
%
%       % multiply vandermonde matrix and coeffs vector
%       coeffs = coeffs(:);
%       B = vdm*coeffs;
%
  end;

  return;


%**************************************************************************%
function P = mris_gradient_nonlin__siemens_legendre(n,X)
%MRIS_GRADIENT_NONLIN__SIEMENS_LEGENDRE  normalize Legendre polynomial with Siemens's convention

  P = legendre(n,X);

  for m=1:n,
    normfact = (-1)^m*sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)));
    P(m+1,:) = normfact*P(m+1,:);
  end;

  return;


  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__spharm_evaluate.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
