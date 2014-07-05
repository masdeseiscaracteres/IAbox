function [U, V, IL]=AlternatingMinLeakage(H,D,options)
% Interference Alignment solution via alternate minimization
%
% Inputs:
% H: Channel matrix indexing coefficients as {rx user, tx user}(rx ant, tx ant)
%
% Outputs:
% U and V.
%
% Reference:
% K. S. Gomadam, V. R. Cadambe, and S. A. Jafar, "A distributed numerical 
% approach to interference alignment and applications to wireless 
% interference networks,” IEEE Trans. Inf. Theory, vol. 57, no. 6, 
% pp. 3309–3322, 2011.


opts.Tol=1e-10; %Mininum interference leakage to exit
opts.MaxIter=5000; %Maximum number of alternating minimization iterations
opts.Verbose=false;

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

K=length(H);

nT=cellfun('size',H(1,:),2)';
nR=cellfun('size',H(:,1),1);
d=diag(D);

if sum(d<=min(nT,nR))<K
    error('Number of antennas at both sides ot the link must be larger than or equal to the number of streams')
end

V=cell(K,1);
U=cell(K,1);

%% Alternate minimization algorithm starts here
IL=nan(1,opts.MaxIter);

if exist('options','var') && isfield(options,'StartingPoint')
    V=options.StartingPoint.V;
    U=options.StartingPoint.U;
    IL(1)=sum(sum(InterferenceLeakage(H,U,V,D)));
else
    for tx=1:K
        %Start with arbitrary precoders
        V{tx}=orth(randn(nT(tx),d(tx))+1j*randn(nT(tx),d(tx))); %V'*V=I_d
    end
end

%% Algorithm starts here

for iter=1:opts.MaxIter
    for rx=1:K
        %Interference covariance matrix
        Q=0;
        for tx=[1:rx-1 rx+1:K]
            Q=Q+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';
        end
        
        %Obtain the subspace that contains the least interference -> Decoders
        [A,S,~]=svd(Q);
        U{rx}=A(:,end-d(rx)+1:end);  % smallest eigenvectors -> interference free subspace
    end
   
    ILtemp=0;
    for tx=1:K
        %Interference covariance matrix
        Q=0;
        for rx=[1:tx-1 tx+1:K]
            Q=Q+H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx};
        end
        
        %Obtain the subspace that contains the least interference -> Decoders
        [A,S,~]=svd(Q);
        V{tx}=A(:,end-d(tx)+1:end);  % smallest eigenvectors -> interference free subspace
        ILtemp=ILtemp+trace(S(end-d(tx)+1:end,end-d(tx)+1:end));
    end
    
    IL(iter+1)=ILtemp;
    
    if mod(iter,opts.Verbose)==0
        fprintf('         Iter.: %3d / %3d,  Interference leakage: %.2e\n',iter,opts.MaxIter,IL(iter+1));
        semilogy(IL(1:iter+1),'b');
        drawnow;
    end
    
    if abs(IL(iter+1))<opts.Tol
        IL=IL(1:iter+1);
        break
    end
    
end
