function [U, V, SINR]=AlternatingMaxSINR(H,D,W,options)
% MaxSINR via alternate minimization
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

opts.Tol=1e-6; %End condition
opts.MaxIter=5000; %Maximum number of alternating minimization iterations
opts.Verbose=false;
opts.Orthogonalize=false; %Ortogonalization step suggested in 
% I. Santamaria, O. Gonzalez, R. W. Heath Jr., and S. W. Peters,
% "Maximum Sum-Rate Interference Alignment Algorithms for MIMO Channels,"
% in 2010 IEEE Global Telecommunications Conference GLOBECOM, 2010, pp. 1–6.

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

SINR=cell(K,1);
%% Alternate minimization algorithm starts here
Cost=nan(sum(diag(D)),opts.MaxIter);

if exist('options','var') && isfield(options,'StartingPoint')
    V=options.StartingPoint.V;
else
    for tx=1:K
        %Start with arbitrary precoders
        V{tx}=orth(randn(nT(tx),d(tx))+1j*randn(nT(tx),d(tx))); %V'*V=I_d
    end
end

%% Algorithm starts here

for iter=1:opts.MaxIter
    %% Calculate decoders
    for rx=1:K
        %Interference covariance matrix at rx_st-th stream of the rx-th receiver
        Q=W{rx}; %Noise covariance matrix
        for tx=1:K
            Q=Q+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}'/D(tx,tx);
        end
        
        for rx_st=1:D(rx,rx)
            R=H{rx,rx}*V{rx}(:,rx_st)*V{rx}(:,rx_st)'*H{rx,rx}'/D(rx,rx);
            B=Q-R;
            
            %Obtain solution to the generalized eigenvalue problem
            % Option 1:
            Uaux=B\(H{rx,rx}*V{rx}(:,rx_st));
            U{rx}(:,rx_st)=Uaux/norm(Uaux);
            SINR{rx}(rx_st)=abs(U{rx}(:,rx_st)'*R*U{rx}(:,rx_st)/(U{rx}(:,rx_st)'*B*U{rx}(:,rx_st)));
            % Option 2:
            %[U_2 SINR_2]=eig(R,B);
            %SINR{rx}(rx_st)=max(diag(abs(SINR_2)));
        end
        
        %Optional step to enforce orthogonality
        if opts.Orthogonalize
            [QU,unused]=qr(U{rx});
            U{rx}=QU(:,1:D(rx,rx));
        end
    end
    
    %% Calculate precoders
    for tx=1:K
        %Interference covariance matrix at rx_st-th stream of the rx-th receiver
        Q=W{tx}; %Noise covariance matrix
        for rx=1:K
            Q=Q+H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx}/D(rx,rx);
        end
        
        for tx_st=1:D(tx,tx)
            R=H{tx,tx}'*U{tx}(:,tx_st)*U{tx}(:,tx_st)'*H{tx,tx}/D(tx,tx);
            B=Q-R;
            Vaux=B\(H{tx,tx}'*U{tx}(:,tx_st));
            V{tx}(:,tx_st)=Vaux/norm(Vaux);
        end
        
        %Optional step to enforce orthogonality
        if opts.Orthogonalize
            [QV,unused]=qr(V{tx});
            V{tx}=QV(:,1:D(tx,tx));
        end
    end
    
    %% Evaluate cost function
    Cost(:,iter)=cat(2,SINR{:})';
    if mod(iter,opts.Verbose)==0
        fprintf('         Iter.: %3d / %3d,  Avg. SINR: %.2e\n',iter,opts.MaxIter,mean(Cost(:,iter)));
        plot(Cost(:,1:iter)','b');
        drawnow;
    end
    
    if iter>1
        if norm(Cost(:,iter)-Cost(:,iter-1))<opts.Tol
            SINR=Cost(:,1:iter);
            break
        end
    end
end
SINR=Cost;