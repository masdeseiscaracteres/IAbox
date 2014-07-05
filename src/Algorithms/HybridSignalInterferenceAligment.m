function [U, V, IL]=HybridSignalInterferenceAligment(H,D,options)
% Interference Alignment solution via joint signal and interference
% optimization. Thia algorithm iteratively reduces both the interference
% power that "leaks" into the signal subspace, and the signal power that
% "leaks" into the interference subspace.
%
% Inputs:
% H: Channel matrix indexing coefficients as {rx user, tx user}(rx ant, tx ant)
%
% Outputs:
% U: Decoders
% V: Precoders
% IL: Interference leakage
%
% Reference:
% K. Raj Kumar and F. Xue, "An iterative algorithm for joint signal and
% interference alignment,” in 2010 IEEE International Symposium on
% Information Theory, 2010, pp. 2293–2297.
%
% Óscar González Fernández
% Advanced Signal Processing Group (GTAS)
% Department of Communications Engineering (DICOM)
% University of Cantabria

opts.Tol=1e-10; %Mininum interference leakage to exit
opts.MaxIter=5000; %Maximum number of alternating minimization iterations
opts.Verbose=false;
opts.w=0.002; %w=1: maximum signal, w=0: minimum interference
opts.EndingUpdate=1e-3;

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


%Algorithm starts here

for iter=1:opts.MaxIter
    for rx=1:K
        %Interference covariance matrix
        Q=-opts.w*H{rx,rx}*V{rx}*V{rx}'*H{rx,rx}';
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
        Q=opts.w*H{tx,tx}'*(eye(nR(tx))-U{tx}*U{tx}')*H{tx,tx};
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
    
    if iter>2
        s=0;
        for us=1:K
            s=s+norm(Vprev{us}-V{us}*V{us}'*Vprev{us},'fro');
            s=s+norm(Uprev{us}-U{us}*U{us}'*Uprev{us},'fro');
        end
        if s<opts.EndingUpdate
            IL=IL(1:iter+1);
            break;
        end
    end
    Vprev=V;
    Uprev=U;
    
    if abs(IL(iter+1))<opts.Tol
        IL=IL(1:iter+1);
        break
    end
    
end
