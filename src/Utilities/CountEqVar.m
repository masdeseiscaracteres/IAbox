function [Ne Nv]=CountEqVar(nT, nR, d, K,int_graph)
% This function counts IA equations and free variables in an 
% interference channel

nT=nT(:);
nR=nR(:);
d=d(:);

%If interference graph matrix is not passed as an argument, assume the
%network is fully connected
if ~exist('int_graph','var')
   int_graph=double(~eye(K,K)); 
end

txs=diag(any(int_graph,2)); %Active transmitters
rxs=diag(any(int_graph)); %Active receivers

Nv=d'*txs*nT+d'*rxs*nR-d'*(txs+rxs)*d;%Count variables
Ne=d'*int_graph*d;%Count equations

%% %ToDo: Generalize to partially connected X networks
% L=length(nT); %Number of transmitters
% K=length(nR); %Number of receivers

% D=ones(K,L); %Demands matrix
% A=ones(K,L); %Adjacency matrix
% cT=sum(D,1).'; %Number of columns
% cR=sum(D,2); %Number of rows


%Double check with the correct values for a symmetric scenario
%Nv=L*K*d(1,1)*(nT(1)+nR(1)-(L+1)*d(1,1))
%Ne=L^2*K*d(1,1)^2*(K-1)

return
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