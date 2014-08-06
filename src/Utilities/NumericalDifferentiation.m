function [DF_X DF_Xconj]=NumericalDifferentiation(fun,X,delta)
% Numerical complex differentiation of scalar complex functions
if nargin==2
delta=1e-4;
end
sz=size(X);
DF_X=nan(sz);
DF_Xconj=nan(sz);
for kk=1:numel(X)
    DX=zeros(sz);
    DX(kk)=delta;
    DR=(fun(X+DX)-fun(X))/delta;
    DX=zeros(sz);
    DX(kk)=1j*delta;
    DI=(fun(X+DX)-fun(X))/delta;
    DF_X(kk)= 0.5*(DR-1j*DI);
    DF_Xconj(kk)= 0.5*(DR+1j*DI);
end