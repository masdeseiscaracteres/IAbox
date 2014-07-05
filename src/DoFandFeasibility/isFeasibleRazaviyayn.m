function feasible=isFeasibleRazaviyayn(nT,nR,D,options)
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
% M. Razaviyayn, G. Lyubeznik, and Z.-Q. Luo, "On the degrees of freedom 
% achievable through interference alignment in a MIMO interference channel,"
% IEEE Transactions on Signal Processing, vol. 60, no. 2, pp. 812–821, Feb. 2012.

%% Default options
opts.Adj=eye(length(nT)); %Adjacency matrix of a fully connected scenario

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%% Prepare input variables
nT=nT(:);
nR=nR(:);
if isvector(D)
    d=D(:);
else
    d=diag(D);
end
K=length(nT);
%%
feasible=-1; %Default output: we know nothing

%% Cases when necessary and sufficient conditions align
if all(d(1)==d) % Every user transmits the same number of streams
    if all(nT(1)==nT) && all(nR(1)==nR) &&  p2pBounds(nT,nR,d,opts)~=0 && ((mod(nT(1),d(1))==0) || (mod(nR(1),d(1))==0))
        feasible=(nT(1)+nR(1))>=(K+1)*d(1);
        return;
    elseif all(mod(nT,d(1))==0) && all(mod(nR,d(1))==0) && (allSubsetsBounds(nT,nR,d,opts)==-1)
        feasible=1;
        return;
    end
end

%In this case only necessary conditions are provided
if p2pBounds(nT,nR,d,opts)==0 || pairwiseBounds(nT,nR,d,opts)==0 || allSubsetsBounds(nT,nR,d,opts)==0
    feasible=0;
end
end

%% O(K) P2P conditions
function feasible=p2pBounds(nT,nR,d,opts)
feasible=-1;
if any(min(nT,nR)<d);
    feasible=0;
end
end

%% O(Neq) conditions
function feasible=pairwiseBounds(nT,nR,d,opts)
feasible=-1;
[rxs, txs, values]=find(opts.Adj); %List of interference links
if any(max(nT(txs),nR(rxs))<(d(txs)+d(rxs)))
    feasible=0;
end
end

%% O(2^Neq) conditions
function feasible=allSubsetsBounds(nT,nR,d,opts)
feasible=-1;

varsT=(nT-d).*d;
varsR=(nR-d).*d;

eqs=(d*d').*not(opts.Adj);

[rxs, txs, values]=find(opts.Adj); %List of interference links
L=length(values);

for subset_idx=1:2^L %All non-empty subsets (2^Neq-1)
    subset_sel=logical(dec2bin(subset_idx,L)-'0');
    tx_sel=unique(txs(subset_sel));
    rx_sel=unique(rxs(subset_sel));
    
    Nv=sum(varsT(tx_sel))+sum(varsR(rx_sel));
    Ne=sum(sum(eqs(rx_sel,tx_sel)));
    if Ne>Nv
        feasible=0;
        return;
    end
end
end