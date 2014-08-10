function [Ne, Nv]=CountEqVar(nT,nR,D,opts)
% COUNTEQVAR Counts the equations and free variables involved in the 
% polynomial equations describing the IA problem for the system described 
% by the input arguments nT, nR and D, where nT (nR) is a vector of
% transmit (receive) antennas and D is the demands matrix, i.e., D(i,j) 
% denotes the number of streams the j-th transmitter wishes to send to the
% i-th receiver.
%
% Usage:
%
% [Ne, Nv]=CountEqVarX(nT,nR,D)
% 
% This is a generalization of the results in the reference below to 
% partially connected networks. The following assumptions have been made:
% - A zero-gain channel does not impose any constraint (or equation).
% - A node which is either not causing interference or being interfered is
%   ignored when counting free variables.
%
% Reference: 
% Sun, Hua, Tiangao Gou, and Syed Ali Jafar "Degrees of freedom of MIMO X 
% networks: Spatial scale invariance and one-sided decomposability,"
% IEEE Transactions on Information Theory, vol. 59, no. 12, 8377-8385.

nT=nT(:); %Number of TX antennas as a column vector
nR=nR(:); %Number of RX antennas as a column vector
[N, M] = size(D); %Number of receivers and transmitters

%% Number of variables
sumDtx=sum(D,2);
Nv_tx=sum(sum(bsxfun(@minus,nT',D).*D));
Nv_rx=(nR-sumDtx)'*sumDtx;
Nv=Nv_tx+Nv_rx;

%% Number of equations
Ne=0;
for rx=1:N
Ne=Ne+sum(sumDtx([1:rx-1 rx+1:N]))*sum(D(rx,:));
end

%It may be useful...

%WRONG! Nv=cT.'*(nT-cT)+(nR-cR).'*cR; %Number of variables
%WRONG! Ne=sum(sum(A.*(cR*cT.'))); %Number of equations
% [lmat kmat]=meshgrid(1:L,1:K);
% 
% Dp=arrayfun(@(k,l) sum(D([1:k-1 k+1:end],l)),kmat,lmat);
% Arow_part=(A.*Dp)~=0; %Arow_part is 1 if there are equations involving Hkl
% 
% Ap=arrayfun(@(k,l) sum(A([1:k-1 k+1:end],l)~=0),kmat,lmat);
% Acol_part=and(D,Ap); %Acol_part(j,p) is 1 if there are equations involving Vjp
% 
% Act_rx = find(sum(A.*D,2)>0); % Active_Rx_finder