function Qdec2 = rmQdecCol(Qdec1,col)
% Qdec2 = rmQdecCol(Qdec1,col)
%
% Removes the specify column from cell string array Qdec1. 
%
% Input
% Qdec1: Two dimensional cell string array of Qdec data (eg. read with 
% fReadQdec).
% col: Column to be removed.
%
% Output
% Qdec2: Two dimensional cell string array of Qdec data.
%
% $Revision: 1.2 $  $Date: 2012/12/12 22:58:13 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: vinke $
%    $Date: 2012/12/12 22:58:13 $
%    $Revision: 1.2 $
%
if nargin < 2
    error('Too few inputs');
end;
Qdec2 = [Qdec1(:,1:col-1) Qdec1(:,col+1:end)]; 


