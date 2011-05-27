function [X, Nc] = fast_par2Xtask(par,ntrs,TR,ERM,W)
% [X, Nc] = fast_par2Xtask(par,ntrs,TR,ERM,W)
%
% ERM - vector of event-response model structures. The
% structure has (at least) the following fields:
%   name - string with the name of the model
%   params - vector with parameter list
%
% If the length of ERM is 1, then the same ERM is applied
% to all event types.
%
% Possible ERMs and their parameters:
%  fir         - tprestim, ter, timewindow
%  gamma       - delay, dispersion, boxcarwidth
%  gamma+deriv - delay, dispersion, boxcarwidth
%
% See also:
%   fast_sched2Xerm, fast_sched2Xfir, 
%   fast_sched2Xgamma, fast_sched2Xgammaderiv, fast_par2nconds


%
% fast_par2Xtask.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:04 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

X = [];
Nc = [];

if(nargin ~= 5)
  fprintf('[X, Nc] = fast_par2Xtask(par,ntrs,TR,ERM,W)\n');
  return;
end

% Get the number of conditions %
[Nc CondList Holes] = fast_par2nconds(par);
if(Holes)
  fprintf('ERROR: holes in the paradigm\n');
  return;
end

if(length(ERM) ~= Nc & length(ERM) ~= 1)
  fprintf('ERROR: Number of ERMs does not equal number of conditions\n');
  return;
end

X0 = [];
for c = 1: Nc
  indc = find(par(:,2)==c);

  tPres = par(indc,1);

  if(~isempty(W)) Wc = W(indc);
  else Wc = [];
  end

  if(length(ERM) ~= 1) ERMc = ERM(c);
  else                 ERMc = ERM;
  end

  Xc = fast_sched2Xerm(tPres,ntrs,TR,Wc,ERMc);
  if(isempty(Xc)) return; end

  X0 = [X0 Xc];
end

X = X0;

return;
