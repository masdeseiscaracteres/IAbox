function [Feasible resid] = isFeasibleNumericalTest(nT,nR,D,options)
%isFeasibleNumericalTest   Interference alignment feasibility test
%   Feasible = isFeasibleNumericalTest(nT,nR,D) checks the feasibility of
%   the interference channel passed in the parameters. Vectors nT and nR
%   contain the number of antennas at the transmitting and receiving side,
%   respectively. The K x K (where K is the number of users) matrix D 
%   contains the number of streams each user wishes to send to its 
%   corresponding receiver in the main diagonal (a vector is also accepted).
%   The length of nT and nR must be equal to K.
%   Returns Feasible = 1 if the system is feasible and Feasible = 0 if not.
%
%   Feasible = isFeasibleNumericalTest(nT,nR,D,struct('A',A)) check the feasibility of
%   a system whose connectivity is specified by the K x K matrix A. If
%   the i-th receiver cannot be reached from the j-th transmitter then
%   A(i,j)=0.
%
%   [Feasible,resid] = isFeasibleNumericalTest(...) also returns a residual
%   indicating the degree of confidence in the result. A large residual
%   means the system is likely to be infeasible whereas a small residual
%   indicates the system is feasible.
%
%   isFeasibleNumericalTest(...,options) uses the structure "options" to pass
%   optional arguments to the function. The only allowed option so far
%   controls the verbosity of the function:
%       options.Verbose = 0: Silent mode
%       options.Verbose = 1: Verbose mode
%       options.Verbose = 2: Provides numerical details
%
%   Usage examples:
%
%   Check the feasibility of the (2x2,1)^3 system:
%    	Feasible = isFeasibleNumericalTest([2 2 2],[2 2 2],diag([1 1 1]));
%
%   Check the feasibility of a random system:
%       MAX_ANT = 10; %Maximum number of antennas per link side
%       MAX_USERS = 6; %Maximum number of users
%       MAX_STREAMS = 3; %Maximum number of streams
%       K = randi([3 MAX_USERS],1); %Actual number of users
%       nT = randi([2 MAX_ANT],1,K); %Number of TX antennas
%       nR = randi([2 MAX_ANT],1,K); %Number of RX antennas
%       D = diag(randi([1 MAX_STREAMS],1,K)); %Number of streams
%       options.A = logical(sprand(K,K,0.8)); %Random connectivity matrix
%       options.Verbose = 1;
%       Feasible = isFeasibleNumericalTest(nT,nR,D,options);
%
%   Reference:
%
%   Ó. González, C. Beltrán and I. Santamaría, "A Feasibility Test for
%   Linear Interference Alignment in MIMO Channels with Constant
%   Coefficients", IEEE Transactions on Information Theory, vol. 60, no. 3,
%   pp. 1840-1856, Mar. 2014.

%% Turn vectors into column vectors
nT=nT(:);
nR=nR(:);
if isvector(D)
    d=D(:);
else
    d=diag(D);
end
K=length(d);

%% 
opts.Verbose=0;
opts.A=ones(K);

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

% List of interfering transmitters and interfered receivers (they may be repeated)
Adj=(~eye(K)).*logical(opts.A);
[phi_R_rep, phi_T_rep]=find(Adj);
phi_R=unique(phi_R_rep);
phi_T=unique(phi_T_rep);

%% Check two assumptions
% Assumption 1: Point-to-point bounds for direct channels have to be
% satisfied
a1=all(d<=nT & d<=nR); %1: assumption is satisfied
% Assumption 2: Each stream has to be cancelled by either the transmitter
% or the receiver
a2=all(nT(phi_T_rep).*nR(phi_R_rep)>d(phi_T_rep).*d(phi_R_rep)); %1: assumption is satisfied

if not(all([a1 a2]))
    if opts.Verbose
        disp('Violating assumptions -> NOT FEASIBLE')
    end
    Feasible=0;
    resid=NaN;
    return;
end

%% Check if the overall number of equations is larger than the number of variables
Ne=sum(sum((d*d.').*Adj));
Nv=(nR(phi_R)-d(phi_R))'*d(phi_R)+(nT(phi_T)-d(phi_T))'*d(phi_T);

if Ne>Nv
    if opts.Verbose
        disp('Improper system -> NOT FEASIBLE')
    end
    Feasible=0;
    resid=NaN;
    return;
end

%% Check if mapping is surjective

%Particular channel representatives
[H, A, B] = randomHAB_I(nT,nR,d,K);
%Build the linear mapping
theta=BuildSparseIAJacobianMapping(nT,nR,d,K,A,B,Adj);

%Solve Nsys systems of linear equations to verify if matrix theta is singular
Nsys=5;
w=randn(Ne,Nsys)+1j*randn(Ne,Nsys);
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:nearlySingularMatrix');

%Enable solver verbosity
if opts.Verbose==2
    spparms('spumoni',2);
else
    spparms('spumoni',0);
end

% Solve Nsys linear systems
x=theta\w;

%Show solver configuration
if opts.Verbose==2
    spparms
end

%Compute residuals
e=sum(abs(theta*x-w).^2);
resid=sqrt(sum(e)); %Frobenius norm

if sum(e<1e-3)>Nsys/2
    Feasible=1;
else
    Feasible=0;
end

if opts.Verbose
    feas_str={'NOT FEASIBLE','FEASIBLE'};
    fprintf('This system is %s\n',feas_str{Feasible+1});
end
warning('on','MATLAB:rankDeficientMatrix');
warning('on','MATLAB:nearlySingularMatrix');
end

function [H A B] = randomHAB_I(nT,nR,d,K)

H=cell(K,K);
A=cell(K,K);
B=cell(K,K);

for tx=1:K
    for rx=1:K
        Z1=zeros(d(rx),d(tx));
        Z2=zeros(nR(rx)-d(rx),nT(tx)-d(tx));
        A{rx,tx}=randn(d(rx),nT(tx)-d(tx))+1j*randn(d(rx),nT(tx)-d(tx));
        B{rx,tx}=randn(nR(rx)-d(rx),d(tx))+1j*randn(nR(rx)-d(rx),d(tx));
        H{rx,tx}=[Z1 A{rx,tx}; B{rx,tx} Z2];
    end
end
end

function theta=BuildSparseIAJacobianMapping(nT,nR,d,K,A,B,Adj)

%Turn vectors into column vectors
nT=nT(:);
nR=nR(:);
d=d(:);


[rxs txs ~]=find(Adj); %Set of tuples denoting interfering links (rx,tx)

aux=(d*d.').*Adj;
%Dimension of the output vector
Ne=sum(aux(:));
%Dimension of the input vector
Nv=sum((nR+nT-2*d).*d);
%Number of non-zero elements (nnz)
% sum(Nk-dk)dldk when link l-k is present + sum(Ml-dl)dkdl when link
% k-l is present in the network
nzmax=sum((nR-d)'*aux+(nT-d)'*aux');

%Sparseness factor or density
% F=nzmax/Ne/Nv;
% fprintf('Sparseness factor: %.3f\n',F);

cols_HV=nR-d;
cols_UtH=nT-d;

cols_BV=d.*cols_HV;
cols_BU=d.*cols_UtH;
offset_side=sum(cols_BU(1:K));

iijj=zeros(nzmax,3);
c=1;

row=1;
for kk=1:length(rxs)
    tx=txs(kk);
    rx=rxs(kk);
    for tx_st=1:d(tx)
        for rx_st=1:d(rx)
            offset_tx=sum(cols_BU(1:tx-1)); %TX refers to the index of Vdot
            offset_rx=sum(cols_BV(1:rx-1))+offset_side; %RX refers to the index of Udot
            offset_tx_st=(tx_st-1)*cols_UtH(tx); %Refers to the column tx_st of Vdot_tx
            offset_rx_st=rx_st-1; %Refers to the column rx_st of Udot_rx
            
            idxs1=offset_tx+offset_tx_st+(1:cols_UtH(tx));
            idxs2=offset_rx+offset_rx_st+(1:d(rx):cols_BV(rx));
            
            %Store indexes and values
            len=cols_UtH(tx)+cols_HV(rx);
            iijj(c:c+len-1,:)=[row*ones(len,1) [idxs1 idxs2]',[A{rx,tx}(rx_st,:)'; B{rx,tx}(:,tx_st)]];
            c=c+len;
            
            %Move on to the next row
            row=row+1;
        end
    end
end

%Build sparse matrix
theta=sparse(iijj(:,1),iijj(:,2),iijj(:,3),Ne,Nv);

end
