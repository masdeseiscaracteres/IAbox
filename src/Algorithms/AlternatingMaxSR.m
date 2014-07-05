function [U, V, Power]=AlternatingMinimizationMaxSR(H,D,var,options)

%%---------------  Max Sum-Rate Interference alignment  --------------%
%
%
% We enforce a Max Sum-Rate IA solution by modifying the alternating
% minimization Jafar's algorithm.
% At each iteration, for each user, we compute the minor subspace of the interference
% covariance matrix. In the new step we move this solution along the
% geodesic in the direction that maximizes the sum-rate.
% The step size is decreased (annealed)
% in such a way that as the iterations proceed we finally obtain
% a perfect IA solution.
%
% Reference:
% I. Santamaría, Ó. González, R. W. Heath Jr., and S. W. Peters,
% "Maximum Sum-Rate Interference Alignment Algorithms for MIMO Channels,"
% in 2010 IEEE Global Telecommunications Conference (GLOBECOM), 2010.

opts.Tol=1e-10; %Mininum interference leakage to exit
opts.MaxIter=5000; %Maximum number of alternating minimization iterations
opts.Tstart=1; %Annealing factor to start with
opts.eta=0.995; %Annealing factor
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

%%---- Initialization (precoders)-------
U=cell(K,1);
V=cell(K,1);

if exist('options','var') && isfield(options,'StartingPoint')
    U=options.StartingPoint.U;
    V=options.StartingPoint.V;
else
    for us=1:K
        %Start with arbitrary precoders and decoders
        U{us}=orth(randn(nR(us),d(us))+1j*randn(nR(us),d(us))); %U'*U=I_d
        V{us}=orth(randn(nT(us),d(us))+1j*randn(nT(us),d(us))); %V'*V=I_d
    end
end

%%
Power=nan(K,K,opts.MaxIter);     % stores the receive power
P=nan(1,opts.MaxIter);
IL=nan(1,opts.MaxIter);

Qrx=cell(K,K); %Interference covariance matrices at each receiver from each transmitter
Qtx=cell(K,K); %Interference covariance matrices at each transmitter from each receiver
%% Alternating minimization starts here
t=opts.Tstart;
for n=1:opts.MaxIter
    t=t*opts.eta;   % We decrease the step size taken along the Grassman manifold
    
    %Compute covariance matrices
    for rx=1:K
        for tx=1:K
            Qrx{rx,tx}=H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';
            Qtx{rx,tx}=H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx};
        end
    end
    
    %% ---- Decoders update ---- &
    for rx=1:K
        Q=sum(cat(3,Qrx{rx,[1:rx-1 rx+1:K]}),3);
        [A,S,B]=svd(Q);
        U{rx}=A(:,nR(rx)-d(rx)+1:end);  % smallest eigenvectors -> least-interf subspace -> Projector U{rx}*U{rx}'
        
        %% ---- Now we move the interference free solution towards the Max-SR direction--%
        
        % Update covariance matrices involving U{rx}
        for tx=1:K
            Qtx{rx,tx}=H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx};
        end
        
        DeltaSR=0;  % derivative with respect to decoder
        for tx=1:K
            Cbar=sum(cat(3,Qtx{[1:tx-1 tx+1:K],tx}),3)+var*eye(nT(tx));
            C=Cbar+Qtx{tx,tx};
            invCbar=inv(Cbar);
            invC=inv(C);
            
            
            pos=H{rx,tx}*invC*H{rx,tx}'-trace(U{rx}'*H{rx,tx}*invC*H{rx,tx}'*U{rx}); %Positive addend
            %Negative addend
            if tx==rx
                neg=0;
            else
                neg=H{rx,tx}*invCbar*H{rx,tx}'+trace(U{rx}'*H{rx,tx}*invCbar*H{rx,tx}'*U{rx});
            end
            DeltaSR=DeltaSR+pos-neg;
        end
        DeltaSR=DeltaSR*U{rx};
        
        %------------------ solution along the geodesic in the Grassman manifold
        DeltaU=(eye(nR(rx))-U{rx}*U{rx}')*DeltaSR;  % gradient of the distance at U in the tangent space
        [Uaux,Eaux,Vaux]=svd(DeltaU,'econ');
        U{rx}=U{rx}*Vaux*diag(cos(diag(Eaux)*t))*Vaux'+Uaux*diag(sin(diag(Eaux)*t))*Vaux';
    end
    
    %% Evaluate cost function
    for rx=1:K
        for tx=1:K
            %------ Interference covariance matrix produced by the tx-th Tx into the rx-th Rx-----%
            C=U{rx}'*H{rx,tx}*V{tx};
            Power(rx,tx,n)=norm(C,'fro')^2;
            %-----------------------------------------------------------------------------------%
        end
    end
    P(n)=trace(Power(:,:,n));
    IL(n)=sum(sum(Power(:,:,n)))-P(n);
    if mod(n,opts.Verbose)==0
    
        fprintf('         Iter.: %3d / %3d,  Interference leakage: %.2e\n',n,opts.MaxIter,IL(n));
        subplot(2,1,1);
        semilogy(P(1:n),'b');
        drawnow;
        subplot(2,1,2);
        semilogy(IL(1:n),'b');
        drawnow;
    end
    
    if abs(IL(n))<opts.Tol
        Power=Power(:,:,1:n);
        break
    end
    
    %% Compute covariance matrices
    for rx=1:K
        for tx=1:K
            Qrx{rx,tx}=H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';
            Qtx{rx,tx}=H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx};
        end
    end
    %% ---- Precoders update (V) ----%
     for tx=1:K
        Q=sum(cat(3,Qtx{[1:tx-1 tx+1:K],tx}),3);
        [A,S,B]=svd(Q);
        V{tx}=A(:,nT(tx)-d(tx)+1:end);  % smallest eigenvectors -> least-interf subspace -> Projector V{tx}*V{tx}'
        
        %% ---- Now we move the interference free solution towards the Max-SR direction--%
        
        % Update covariance matrices involving U{rx}
        for rx=1:K
            Qrx{rx,tx}=H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';
        end
        
        DeltaSR=0;  % derivative with respect to decoder
        for rx=1:K
            Cbar=sum(cat(3,Qrx{rx,[1:rx-1 rx+1:K]}),3)+var*eye(nR(rx));
            C=Cbar+Qrx{rx,rx};
            invCbar=inv(Cbar);
            invC=inv(C);
       
            pos=H{rx,tx}'*invC*H{rx,tx}-trace(V{tx}'*H{rx,tx}'*invC*H{rx,tx}*V{tx}); %Positive addend
            %Negative addend
            if tx==rx
                neg=0;
            else
                neg=H{rx,tx}'*invCbar*H{rx,tx}+trace(V{tx}'*H{rx,tx}'*invCbar*H{rx,tx}*V{tx});
            end
            DeltaSR=DeltaSR+pos-neg;
        end
        DeltaSR=DeltaSR*V{tx};
        
        %------------------ solution along the geodesic in the Grassman manifold
        DeltaV=(eye(nT(tx))-V{tx}*V{tx}')*DeltaSR;  % gradient of the distance at U in the tangent space
        [Uaux,Eaux,Vaux]=svd(DeltaV,'econ');
        V{tx}=V{tx}*Vaux*diag(cos(diag(Eaux)*t))*Vaux'+Uaux*diag(sin(diag(Eaux)*t))*Vaux';
     end 
end
