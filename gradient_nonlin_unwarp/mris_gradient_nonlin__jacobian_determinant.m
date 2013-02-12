function JacDet = mris_gradient_nonlin__jacobian_determinant(Dx, Dy, Dz, d, X, Y, Z)
%MRIS_GRADIENT_NONLIN__JACOBIAN_DETERMINANT
%
% JacDet = mris_gradient_nonlin__jacobian_determinant(Dx, Dy, Dz, d, X, Y, Z)
%
% this returns the jacobian determinant |dF| of the distortion mapping 
%
%     F: T -> M
%
% where "T" is the true object and "M" is the measured object. this is to be
% multiplied by I(M) when estimating I(T).
%
% if mapping is identity, |dF| will be ones everywhere.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2010/oct/09
% $Id: mris_gradient_nonlin__jacobian_determinant.m,v 1.2 2012/03/11 01:16:14 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  %% TODO: regularize Jacobian determinant, check for extreme values, and
  %% consider Siemens's approach (which appears to take derivatives over a
  %% wider area)

  %% TODO: take coordinate orientation into account (e.g., RAS vs LAI)

  if ( [prod(size(X)) == 1] || [prod(size(Y)) == 1] || [prod(size(Z)) == 1] ),
    JacDet = 1;
    warning(sprintf('==> [%s]: scalar data found; skipping jacobian calculation', mfilename));
    return;
  end;
  
  d1 = sqrt( [X(2,1,1)-X(1,1,1)]^2 + [Y(2,1,1)-Y(1,1,1)]^2 + [Z(2,1,1)-Z(1,1,1)]^2 );
  d2 = sqrt( [X(1,2,1)-X(1,1,1)]^2 + [Y(1,2,1)-Y(1,1,1)]^2 + [Z(1,2,1)-Z(1,1,1)]^2 );
  d3 = sqrt( [X(1,1,2)-X(1,1,1)]^2 + [Y(1,1,2)-Y(1,1,1)]^2 + [Z(1,1,2)-Z(1,1,1)]^2 );

%  d1 = d(1);
%  d2 = d(2);
%  d3 = d(3);

  if ( d1 == 0 || d2 == 0 || d3 == 0 ),
    beep
    disp('weirdness found');
    keyboard
  end;

  %    [1 1 2 3 4] / 2
  %  + [1 2 3 4 4] / 2
  % ------------------

  % dDxdx = diff(Dx,1,2)/dx;    dDxdx = (cat(2,dDxdx(:,1,:),dDxdx) + cat(2,dDxdx,dDxdx(:,end,:))) / 2;
  % dDydy = diff(Dy,1,1)/dy;    dDydy = (cat(1,dDydy(1,:,:),dDydy) + cat(1,dDydy,dDydy(end,:,:))) / 2;
  % dDzdz = diff(Dz,1,3)/dz;    dDzdz = (cat(3,dDzdz(:,:,1),dDzdz) + cat(3,dDzdz,dDzdz(:,:,end))) / 2;

  % this is really gradient with respect to R,C,S, not X,Y,Z
  [dDxdx, dDxdy, dDxdz] = gradient(Dx, d1, d2, d3);
  [dDydx, dDydy, dDydz] = gradient(Dy, d1, d2, d3);
  [dDzdx, dDzdy, dDzdz] = gradient(Dz, d1, d2, d3);

  % recall that siemens does not store transformation F but rather stores
  % the deviation from linearity, so that the components of F are given by
  %
  %    Fx = Dx + x
  %    Fy = Dy + y
  %    Fz = Dz + z
  %
  % this means that:
  %
  %    dFx/dx = dDx/dx + 1
  %    dFx/dy = dDx/dy    
  %    dFx/dz = dDx/dz    
  %    dFy/dx = dDy/dx
  %    dFy/dy = dDy/dy + 1    
  %    dFy/dz = dDy/dz    
  %    dFz/dx = dDz/dx
  %    dFz/dy = dDz/dy    
  %    dFz/dz = dDz/dz + 1    
  
  % old "grad_unwarp" method did not include cross-terms:
  %  JacDet = (1+dFxdx).*(1+dFydy).*(1+dFzdz);
  %  JacDet = abs(JacDet);

  % include cross-terms, which are non-zero in general
  JacDet = + [(1+dDxdx) .* (1+dDydy) .* (1+dDzdz)] ...
           - [(1+dDxdx) .*    dDydz  .*    dDzdy]  ...
           - [   dDxdy  .*    dDydx  .* (1+dDzdz)] ...
           + [   dDxdy  .*    dDydz  .*    dDzdx]  ...
           + [   dDxdz  .*    dDydx  .*    dDzdy]  ...
           - [   dDxdz  .* (1+dDydy) .*    dDzdx];

  JacDet =     abs(JacDet);

  MAX_DETERMINANT = 10;
  JacDet( find(JacDet > MAX_DETERMINANT) ) = MAX_DETERMINANT;


  return;


  %  syms dDxdx dDydy dDzdz dDxdy dDydz dDzdx dDxdz dDydx dDzdy
  %  D = [1+dDxdx dDxdy dDxdz; dDydx 1+dDydy dDydz; dDzdx dDzdy 1+dDzdz]
  %  det(D)



  %************************************************************************%
  %%% $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/mris_gradient_nonlin__jacobian_determinant.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
