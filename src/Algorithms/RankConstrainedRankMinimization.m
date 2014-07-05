function RankConstrainedRankMinimization(H,D,epsilon,MaxIter,IleakTol)

% Reference:
% D. S. Papailiopoulos and A. G. Dimakis, "Interference Alignment as a 
% Rank Constrained Rank Minimization," IEEE Trans. Signal Process., vol. 60,
% no. 8, pp. 4278–4288, Aug. 2012.

%% Parameters
cvx_quiet(true);
[~, K]=size(H);
nT=cellfun('size',H(1,:),2)';
nR=cellfun('size',H(:,1),1);

%% Initialization
V=cell(K,1);
for tx=1:K
    V{tx}=crandn(nT(tx),D(tx,tx));
end

%% Alternating minimization starts here
Ileak=1;
niter=0;
while Ileak>IleakTol
    
    %% Decoder update
    cvx_begin;
        variable Ua(max(nR),max(max((D))),K) complex;
        CostFunction=cvx(0);
        J=cvx([]);
        for rx=1:K
            J=[];
            for tx=1:K
                if tx~=rx
                    J=[J Ua(1:nR(rx),1:D(rx,rx),rx)'*H{rx,tx}*V{tx}]; % Concatenation of all interference channels seen by receiver rx
                end
            end
            CostFunction=CostFunction+norm_nuc(J);
        end
        minimize( CostFunction )
        subject to
        for user=1:K
            real(lambda_min(Ua(1:nR(rx),1:D(user,user),user)'*H{user,user}*V{user})) >= epsilon; % Minimum eigenvalue constraint
            Ua(1:nR(rx),1:D(user,user),user)'*H{user,user}*V{user} == semidefinite(D(user,user)); 
        end
    cvx_end
    for user=1:K
        U{user}=Ua(1:1:nR(user),1:D(user,user),user);
    end        
    
    %% Precoder update    
    cvx_begin;
        variable Va(max(nT),max(max((D))),K) complex;
        CostFunction=cvx(0);
        J=cvx([]);
        for rx=1:K
            J=[];
            for tx=1:K
                if tx~=rx
                    J=[J U{rx}'*H{rx,tx}*Va(1:nT(tx),1:D(tx,tx),tx)]; % Concatenation of all interference channels seen by receiver rx
                end
            end
            CostFunction=CostFunction+norm_nuc(J);
        end
        minimize( CostFunction )
        subject to
        for user=1:K
            real(lambda_min(U{user}'*H{user,user}*Va(1:nT(tx),1:D(user,user),user))) >= epsilon; % Minimum eigenvalue constraint
            U{user}'*H{user,user}*Va(1:nT(tx),1:D(user,user),user) == semidefinite(D(user,user)); 
        end
    cvx_end
    for user=1:K
        V{user}=Va(1:nT(user),1:D(user,user),user);
    end   
    
    %% Interference leakage
    Vorth=V;
    Uorth=U;
    for k=1:K
        [Q,R]=qr(V{k});
        Vorth{k}=Q(:,1:D(k,k));
        [Q,R]=qr(U{k});
        Uorth{k}=Q(:,1:D(k,k));
    end
    Ileak=0;
    for rx=1:K
         for tx=[1:rx-1 rx+1:K]
            Ileak=Ileak+norm(Uorth{rx}'*H{rx,tx}*Vorth{tx},'fro');
         end
    end
    Ileak
    niter=niter+1;
    if niter==MaxIter
        break;
    end

end
V=Vorth;
U=Uorth;

