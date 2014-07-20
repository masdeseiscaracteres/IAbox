function sumD=MaxDoF2UserXChannel(nT,nR)
%MaxDoF2UserXChannel Upper bound on the number of information-theoretic DoF
%   in a 2-user X channel. This bound is proven to be tight if all users
%   are equipped with the same number of antennas.
%
% Reference:
%
% S. A. Jafar and S. Shamai, "Degrees of Freedom Region of the MIMO X Channel,"
% IEEE Transactions on Information Theory, vol. 54, no. 1, Jan. 2008.

bounds=zeros(7,1);
bounds(1)=nT(1)+nT(2);
bounds(2)=nR(1)+nR(2);
bounds(3)=(max([nT(1) nR(1)])+max([nT(1) nR(2)])+nT(2))/2;
bounds(4)=(max([nT(2) nR(1)])+max([nT(2) nR(2)])+nT(1))/2;
bounds(5)=(max([nT(1) nR(1)])+max([nT(2) nR(1)])+nR(2))/2;
bounds(6)=(max([nT(1) nR(2)])+max([nT(2) nR(2)])+nR(1))/2;
bounds(7)=(max([nT(1) nR(1)])+max([nT(1) nR(2)])+max([nT(2) nR(1)])+max([nT(2) nR(2)]))/3;
sumD=min(bounds);