function [U,V]=AlternatingMinLeakageStructuredChannel(H,D,options)

% Reference:
% C. Lameiro, Ó. González, and I. Santamaría, "An Interference Alignment 
% Algorithm for Structured Channels," in IEEE 14th Workshop on Signal 
% Processing Advances in Wireless Communications (SPAWC), 2013, pp. 295--299.


%% Manage options
opts.MaxIter=5000; %Maximum number of iterations
opts.Tol=1e-12; %Target interference leakage
opts.epsilon=1e-3; % Minimum eigenvalue constraint
opts.Improper=false; %Assume proper signaling by default
opts.L=1; %Number of symbols extension (assume no symbol extension by default)
opts.Verbose=false;

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%%
K=length(H);

M=cellfun('size',H(1,:),2)'; %Number of columns in channel matrices
N=cellfun('size',H(:,1),1); %Number of rows in channel matrices
d=diag(D);

%Abbreviate some variable names
L=opts.L;
% nT=M/L; %Obtain number of TX antennas from the dimension of the channel matrices
% nR=N/L; %Obtain number of RX antennas from the dimension of the channel matrices

%% Initialitation

V=cell(K,1);
U=cell(K,1);
for k=1:K
    if opts.Improper
        V{k}=randn(M(k),d(k));   % random real precoder
        U{k}=randn(N(k),d(k));   % random real decoder
    else
        V{k}=crandn(M(k),d(k));   % random complex precoder
        U{k}=crandn(N(k),d(k));   % random complex decoder
    end
    [Q,~]=qr(V{k});
    V{k}=Q(:,1:d(k))/sqrt(d(k));
    [Q,~]=qr(U{k});
    U{k}=Q(:,1:d(k))/sqrt(d(k));
end

%% The algorithm starts here

Cost=zeros(opts.MaxIter+1,1);
Cost(1)=1;
n=2;
while n<=opts.MaxIter+1 && Cost(n-1)>=opts.Tol
        
    %% Decoder update (U)
    for rx=1:K
        Qi=0;
        for tx=1:K
            if tx~=rx
                Qi=Qi+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}'; % Interference covariance matrix
            else
                Qd=H{rx,rx}*V{rx}*V{rx}'*H{rx,rx}'; % Desired signal covariance matrix
            end
        end
        [A,D]=eig(Qi,Qd); % Solve generalized eigenvalue problem
        [eigen,i]=sort(diag(abs(D)),'ascend');
        A=A(:,i);
        Ux=A(:,1:d(rx));
        if opts.Improper & imag(Ux)
            continue;
        end
        Qx=Ux'*Qd*Ux;
        Ux=Ux*Qx^-(1/2)*sqrt(opts.epsilon);
        if trace(Ux'*Ux)>1 % i.e., if constraint is active
            muU=1;
            muL=0;
            %% Bisection method
            while muU-muL>=1e-10
                mu=(muU+muL)/2;
                [A,D]=eig((1-mu)*Qi+mu*eye(length(Qi)),Qd);
                [eigen,i]=sort(diag(abs(D)),'ascend');
                A=A(:,i);
                Ux=A(:,1:d(rx));
                if opts.Improper & imag(Ux)
                    continue;
                end
                Qx=Ux'*Qd*Ux;
                Ux=Ux*Qx^-(1/2)*sqrt(opts.epsilon);
                if trace(Ux'*Ux)>1
                    muL=mu;
                else
                    muU=mu;
                end
            end
        end
        U{rx}=Ux/sqrt(trace(Ux'*Ux));
    end
    
    
    %% Precoder update (V)
    for tx=1:K
        Ri=0;
        for rx=1:K
            if rx~=tx
                Ri=Ri+H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx}; % Interference covariance matrix
            else
                Rd=H{tx,tx}'*U{tx}*U{tx}'*H{tx,tx}; % Desired signal covariance matrix
            end
        end
        [A,D]=eig(Ri,Rd); % Solve generalized eigenvalue problem
        [eigen,i]=sort(diag(abs(D)),'ascend');
        A=A(:,i);
        Vx=A(:,1:d(tx));
        if opts.Improper & imag(Vx)
            continue;
        end
        Rx=Vx'*Rd*Vx;
        Vx=Vx*Rx^-(1/2)*sqrt(opts.epsilon);
        if trace(Vx'*Vx)>1
            muU=1;
            muL=0;
            while muU-muL>=1e-10
                mu=(muU+muL)/2;
                [A,D]=eig((1-mu)*Ri+mu*eye(length(Ri)),Rd);
                [eigen,i]=sort(diag(abs(D)),'ascend');
                A=A(:,i);
                Vx=A(:,1:d(tx));
                if opts.Improper & imag(Vx)
                    continue;
                end
                Rx=Vx'*Rd*Vx;
                Vx=Vx*Rx^-(1/2)*sqrt(opts.epsilon);
                if trace(Vx'*Vx)>1
                    muL=mu;
                else
                    muU=mu;
                end
            end
        end
        V{tx}=Vx/sqrt(trace(Vx'*Vx));
    end
    
    % Cost function update
    for rx=1:K
        Qi=0;
        for tx=1:K
            if tx~=rx
                Qi=Qi+U{rx}'*H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}'*U{rx}; % Interference covariance matrix
            end
        end
        Cost(n)=Cost(n)+abs(trace(Qi));
    end

    if mod(n,opts.Verbose)==0
        fprintf('         Iter.: %3d / %3d,  Interference leakage: %.2e\n',n,opts.MaxIter,Cost(n));
        semilogy(Cost(1:n),'b');
        drawnow;
    end
    
    n=n+1;
end