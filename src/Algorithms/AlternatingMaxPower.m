function [U, V, Power] = AlternatingMaxPower(H,D,options)

%%---------------  Max-Power Interference alignment --------------%
%
%
% We enforce a Max-Power IA solution by modifying the alternating
% minimization algorithm by Jafar.
% At each iteration, for each user, we compute the minor subspace of the interference
% covariance matrix and the major subspace for the desired covariance
% matrix. These two points are the min-interference and max-power
% solutions. In the new step they are combined by moving the
% min-interference solution towards the max-power solution along the
% geodesic between these two points. The step size is decreased (annealed)
% in such a way that as the iterations proceed we remain closer to the
% min-interference solution.
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
    V=options.StartingPoint.V;
else
    for tx=1:K
        %Start with arbitrary precoders
        V{tx}=orth(randn(nT(tx),d(tx))+1j*randn(nT(tx),d(tx))); %V'*V=I_d
    end
end
%%
Power=nan(K,K,opts.MaxIter);     % stores the receive power
P=nan(1,opts.MaxIter);
IL=nan(1,opts.MaxIter);

%% Alternating minimization starts here
t=opts.Tstart;
for n=1:opts.MaxIter
    
    t=t*opts.eta;   % We decrease the step size taken along the Grassman manifold
    
    %% -------------------------- Decoder update (U) ---------------------------------%    
    for rx=1:K
        Q=0; %Int. convariance matrix at the k-th Rx
         for tx=[1:rx-1 rx+1:K]
            Q=Q+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}'; % Interference covariance matrix produced by the l-th Tx into the k-th Rx
         end
        [A,S,B]=svd(Q);
        U{rx}=A(:,nR(rx)-d(rx)+1:end);  % smallest eigenvectors -> least-interf subspace -> Projector U{rx}*U{rx}'
        
        %% ---- Max-power subspace solution at the k-th Rx--%
        [A,S,B]=svd(H{rx,rx}*V{rx}*V{rx}'*H{rx,rx}');
        Umax=A(:,1:d(rx));   % the largest singular vector form the subspace with most power for the desired user
        %-------------------------------------------------------------
        
        %% Now we move the minimum interference solution towards the max power solution
        %             solution along the geodesic in the Grassman manifold
        DeltaU=(eye(nR(rx))-U{rx}*U{rx}')*(Umax*Umax')*U{rx};  % gradient of the distance at U in the tangent space
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
    
    %% -------------------------------  Precoders update (V) --------------------------------%    
    for tx=1:K
         Q=0; %Int. convariance matrix at the k-th Rx
         for rx=[1:tx-1 tx+1:K]
            Q=Q+H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx}; % Interference covariance matrix produced by the rx-th Rx into the tx-th Tx
         end
        [A,S,B]=svd(Q);
        V{tx}=A(:,nT(tx)-d(tx)+1:end);  % smallest eigenvectors -> least-interf subspace -> Projector U{rx}*U{rx}'
        
        %% ---- Max-power subspace solution at the tx-th Tx--%
        [A,S,B]=svd(H{tx,tx}'*U{tx}*U{tx}'*H{tx,tx});
        Vmax=A(:,1:d(tx));   % the largest singular vector form the subspace with most power for the desired user
        %-------------------------------------------------------------
        
        %% Now we move the minimum interference solution towards the max power solution
        %             solution along the geodesic in the Grassman manifold
        DeltaV=(eye(nT(tx))-V{tx}*V{tx}')*(Vmax*Vmax')*V{tx};  % gradient of the distance at V in the tangent space
        [Uaux,Eaux,Vaux]=svd(DeltaV,'econ');
        V{tx}=V{tx}*Vaux*diag(cos(diag(Eaux)*t))*Vaux'+Uaux*diag(sin(diag(Eaux)*t))*Vaux';
    end
end
