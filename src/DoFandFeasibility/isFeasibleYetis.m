function feasible=isFeasibleYetis(nT,nR,D,options)
% Input:
%   nT: Kx1 vector of transmit antennas
%   nR: Kx1 vector of receive antennas
%   D: KxK diagonal matrix with vector of streams as its main diagonal
%
% Output:
%   feasible=1: system is feasible
%   feasible=0: system is infeasible
%   feasible=-1: system feasibility could not be determined
%
% Reference:
% C. M. Yetis, T. Gou, S. A. Jafar, and A. H. Kayran, "On feasibility of 
% interference alignment in MIMO interference networks,"
% IEEE Trans. Signal Process., vol. 58, no. 9, pp. 4771–4782, 2010.

%% Default options
opts.Adj=~eye(length(nT)); %Adjacency matrix of a fully connected scenario

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%%
nT=nT(:);
nR=nR(:);
if isvector(D)
    d=D(:);
else
    d=diag(D);
end
K=length(nT);

%% Default behavior
feasible=-1;

%% Symmetric systems
if all(nT(1)==nT) && all(nR(1)==nR) && all(d(1)==d)
    proper= (nT(1)+nR(1))>=(K+1)*d(1);
    if d(1)==1
        feasible=proper;
    elseif proper==1;
        feasible=-1;
    end
    if proper==0
        feasible=0;
    end
    return;
end

%% Build stream adjacency matrix from link adjacency matrix
A=cell(K,K);
[rxs,txs,values]=find(opts.Adj);
for kk=1:K^2
    [rx tx]=ind2sub([K K],kk);
    A{rx,tx}=opts.Adj(rx,tx)*ones(d(rx),d(tx));
end
A=cell2mat(A);

%% Number of free variables per TX and RX stream
N=sum(d);
varT_s=nT-d; %available dim per TX and stream
varR_s=nR-d; %available dim per RX and stream
varT=zeros(1,N);
varR=zeros(1,N);
ii=1;
for kk=1:K
    for kkk=1:d(kk)
        varT(ii)=varT_s(kk); %available dim per TX
        varR(ii)=varR_s(kk); %available dim per RX
        ii=ii+1;
    end
end
%% O(2^Neq) complexity, exactly 2^Neq-1
[rx_strs,tx_strs]=find(A);
if length(rx_strs)>15
    warning('The number of conditions to evaluate is very large, it may take a while');
end
for strlink_sel_dec=1:2^length(rx_strs)
    strlink_sel=logical(dec2bin(strlink_sel_dec,length(rx_strs))-'0');
    rx_str=unique(rx_strs(strlink_sel));
    tx_str=unique(tx_strs(strlink_sel));
    
    Nv=sum([varR(rx_str) varT(tx_str)]);
    Ne=sum(strlink_sel);
    if Ne>Nv
        %       proper=0; %Improper implies infeasible
        feasible=0;
        return;
    end
end

%%
% If we have reached this point, we do not known anything about the feasibility of the system
% unless every user is transmitting a single stream, in that case the
% system is feasible
if all(d==1)
    feasible=1;
end
