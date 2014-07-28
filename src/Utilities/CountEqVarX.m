function [Ne, Nv]=CountEqVarX(nT,nR,D,opts)

nT=nT(:);
nR=nR(:);

Ne=0;
Nv=0;


% sum(sum(bsxfun(@minus,nT,D).*D))
% 
% sumDtx=sum(D,2);
% (nR-sumDtx)*sumDtx
% 
% keyboard