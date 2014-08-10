function [Ne, Nv]=CountEqVar(nT,nR,D,options)
% COUNTEQVAR Counts the equations and free variables involved in the
% polynomial equations describing the IA problem for the system described
% by the input arguments nT, nR and D, where nT (nR) is a vector of
% transmit (receive) antennas and D is the demands matrix, i.e., D(i,j)
% denotes the number of streams the j-th transmitter wishes to send to the
% i-th receiver.
%
% Usage:
%
% [Ne, Nv]=CountEqVar(nT,nR,D)
%
% This is a generalization of the results in the reference below to
% partially connected networks. The following assumptions have been made:
% - A zero-gain channel does not impose any constraint (or equation).
% - A node which is either not causing interference or being interfered is
%   ignored when counting free variables.
% - Transmitters/receivers are not aware of the network connectivity.
%   Consequently, they cannot adapt the number of columns in their
%   precoders/decoders to that supported by the network. Precoder/decoder
%   dimensions are dictated by the demands matrix.
%
% Reference:
% Sun, Hua, Tiangao Gou, and Syed Ali Jafar "Degrees of freedom of MIMO X
% networks: Spatial scale invariance and one-sided decomposability,"
% IEEE Transactions on Information Theory, vol. 59, no. 12, 8377-8385.

%% Default value for optional parameters
opts.A=ones(size(D)); %Connectivity matrix: fully-connected matrix by default
opts.ConnectivityAware=false; %Transmitters/receivers are not aware of the network connectivity

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%% Some definitions
nT=nT(:); %Number of TX antennas as a column vector
nR=nR(:); %Number of RX antennas as a column vector
[K, L] = size(D); %Number of receivers and transmitters

if opts.ConnectivityAware
    Deff=(opts.A.*D); %Effective demands matrix, if nodes are aware of network connectivity they can adapt their demands
else
    Deff=D;
end
% cT=sum(Deff,1)'; %Number of columns of the precoders
cR=sum(Deff,2); %Number of columns of the decoders

[lmat kmat]=meshgrid(1:L,1:K);
%% Number of variables
% Precoder variables
Ap=arrayfun(@(k,l) sum(opts.A([1:k-1 k+1:end],l)~=0),kmat,lmat);
activeV=and(Deff,Ap); %activeV(j,p) is 1 if there are equations involving Vjp
Nv_tx_jp=bsxfun(@minus,nT',Deff).*Deff; %Every column in a subprecoder is orthogonal to each other but not necessarily to every other subprecoder's columns
Nv_tx=sum(sum(Nv_tx_jp(activeV)));

% Decoder variables
Dp=arrayfun(@(k,l) sum(Deff([1:k-1 k+1:end],l)),kmat,lmat);
temp=sum(opts.A.*Dp,2); %temp(i): Number of interfering streams at the i-th receiver
activeU=(temp>0); % activeU(k) is 1 if there are equations involving Uk
Nv_rx_k=(nR-cR).*cR;%Every column in a decoder is orthogonal to each other
Nv_rx=sum(Nv_rx_k(activeU));

% Total number of variables
Nv=Nv_tx+Nv_rx;

%% Number of equations
Ne=temp'*cR;
