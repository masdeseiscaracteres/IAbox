function feasible=isFeasibleBresler(nT,nR,D,options)
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
% G. Bresler, D. Cartwright, and D. Tse, "Interference alignment for the
% MIMO interference channel," ArXiv preprint available
% http//arxiv.org/abs/1303.5678, Mar. 2013.

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


%%
feasible=-1; %Default output

if K>=3 && all(nT(1)==[nT;nR]) && all(d(1)==d)
    feasible=(2*nT(1)>=(K+1)*d(1));
else
    feasible=allSubsetsBounds(nT,nR,d,opts);
end

end
%% O(2^K) conditions
function feasible=allSubsetsBounds(nT,nR,d,opts)
feasible=-1;
K=length(d);

[rxs, txs, values]=find(opts.Adj); %List of interference links

varsT=(nT-d).*d;
varsR=(nR-d).*d;
eqs=(d*d').*opts.Adj;

for subset_idx=1:2^K
    subset_sel=logical(dec2bin(subset_idx,K)-'0');
    Nv=sum(varsT(subset_sel))+sum(varsR(subset_sel));
    Ne=sum(sum(eqs(subset_sel,subset_sel)));
    if Ne>Nv
        feasible=0;
        return;
    end
end
end
