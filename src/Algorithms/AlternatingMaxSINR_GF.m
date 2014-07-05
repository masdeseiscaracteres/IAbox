function [U, V, SINR]=AlternatingMaxSINR_GF(H,D,W,options)
% MaxSINR via alternate minimization
%
% Inputs:
% H: Channel matrix indexing coefficients as {rx user, tx user}(rx ant, tx ant)
%
% Outputs:
% U and V.
%
% Reference:
% C. M. Yetis, Y. Zeng, K. Anand, Y. L. Guan, and E. Gunawan,
% "Sub-Stream Fairness and Numerical Correctness in MIMO
% Interference Channels," in 3rd IEEE Symposium on Wireless Technology &
% Applications (ISWTA), sep. 2013

opts.Tol=1e-6; %End condition
opts.MaxIter=5000; %Maximum number of alternating minimization iterations
opts.Verbose=false;
opts.Orthogonalizing=false;

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
        Q=0;
        for tx=[1:rx-1 rx+1:K]
            Q=Q+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}'/D(tx,tx); %To-Do: is dividing by D(tx,tx) right? think of assymetric scenarios
        end
        
        R=H{rx,rx}*V{rx}*V{rx}'*H{rx,rx}'/D(rx,rx);
        B=Q+W{rx};
        
        %Obtain solution to the generalized eigenvalue problem
        [U_2 SINR_2]=eig(R,B);
        [values,indices]=sort(diag(abs(SINR_2)),'descend');
        U{rx}=U_2(:,indices(1:D(rx,rx)));
        
        %Optional step to enforce orthogonality
        if opts.Orthogonalize
            [QU,unused]=qr(U{rx});
            U{rx}=QU(:,1:D(rx,rx));
        else
            U{rx}=U{rx}/norm(U{rx},'fro')*sqrt(D(rx,rx));
        end
        
        %Calculate cost function
        for rx_st=1:D(rx,rx)
            SINR{rx}(rx_st)=real((U{rx}(:,rx_st)'*R*U{rx}(:,rx_st))/(U{rx}(:,rx_st)'*B*U{rx}(:,rx_st)));
        end
    end
    
    %% Calculate precoders
    for tx=1:K
        %Interference covariance matrix at rx_st-th stream of the rx-th receiver
        Q=0;
        for rx=[1:tx-1 tx+1:K]
            Q=Q+H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx}/D(rx,rx);
        end
        
        R=H{tx,tx}'*U{tx}*U{tx}'*H{tx,tx}/D(tx,tx);
        B=Q+W{tx};
        
        %Obtain solution to the generalized eigenvalue problem
        [V_2 SINR_2]=eig(R,B);
        [values,indices]=sort(diag(abs(SINR_2)),'descend');
        V{tx}=V_2(:,indices(1:D(tx,tx)));
        
        %Optional step to enforce orthogonality
        if opts.Orthogonalize
            [QV,unused]=qr(V{tx});
            V{tx}=QV(:,1:D(tx,tx));
         else
            V{tx}=V{tx}/norm(V{tx},'fro')*sqrt(D(tx,tx));
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