function [U,V,IL]=SteepestDescentMinLeakage(H,D,options)

opts.MaxIter=1000;
opts.Tol=1e-8;
opts.MinStepSize=1e-5;
opts.Manifold='Grassmann'; %Allowed valueas are: Stiefel, Grassmann or none

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%%
%ToDo: check that D describes an interference channel
[~, K]=size(H);
nT=cellfun('size',H(1,:),2)';
nR=cellfun('size',H(:,1),1);
d=diag(D);

if isfield(opts,'StartingPoint')
    V=opts.StartingPoint.V;
else
    V=arrayfun(@(a,b)orth(crandn(a,b)),nT',d','UniformOutput',false);
end

beta=1*ones(K,1); %Initial step size


IL=nan(1,opts.MaxIter);

IL(1)=CostFunction(H,V,D);
  
switch lower(opts.Manifold)
    case 'stiefel'
        g=@(Z1,Z2,X) real(trace(Z1'*(eye(size(X,1))-0.5*X*X')*Z2)); %Distance metric (inner product)
        egrad2rgrad=@(grad,X) -X*grad'*X+grad; %Euclidean gradient to Riemannian gradient
        retraction=@proj_stiefel; %Retract back to the manifold
    case 'grassmann'
        g=@(Z1,Z2,X) trace(Z1'*Z2);
        egrad2rgrad=@(grad,X) (eye(size(X,1))-X*X')*grad;
        retraction=@proj_grassmann;
    otherwise
        g=@(Z1,Z2,X) trace(Z1'*Z2);
        egrad2rgrad=@(grad,X) grad;
        retraction=@(X) orth(X);
end

for nn=1:opts.MaxIter
    for tx=1:K
        %Compute covariance matrices
        Q=cell(K,1);
        M=cell(K,1);
        for rx=1:K
            Q{rx}=0;
            for tx2=[1:rx-1 rx+1:K]
                Q{rx}=Q{rx}+H{rx,tx2}*V{tx2}*V{tx2}'*H{rx,tx2}';
            end
            
            lambdas=svd(Q{rx});
            M{rx}=0;
            for s=nR(rx)-D(rx,rx)+1:nR(rx)
                M_s=1;
                for kk=[1:s-1 s+1:nR(rx)]
                    M_s=M_s*(lambdas(kk)*eye(nR(rx))-Q{rx})/(lambdas(kk)-lambdas(s));
                end
                M{rx}=M{rx}+M_s;
            end
            
        end
        
        %% Steepest descent direction
        D_V=0;
        for rx=[1:tx-1 tx+1:K]
            D_V=D_V+H{rx,tx}'*M{rx}*H{rx,tx}*V{tx};
        end
        
        %Verify derivative numerically
%         fun=@(X)CostFunction(H,[V(1:tx-1) X V(tx+1:end)],D);
%         [DF_X DF_Xconj]=NumericalDifferentiation(fun,V{tx});
%         assert(norm(D_Vconj-DF_Xconj,'fro')^2<1e-5);
        
        Z=-egrad2rgrad(D_V,V{tx}); %convert from Euclidean gradient to Riemannian gradient
        
        %% Armijo's rule
        B1=V;
        B2=V;
        m=g(Z,Z,V{tx}); %calculate distance on the manifold
        
        if nn>1
            while 1 %Step 3
                B1{tx}=retraction(V{tx}+2*beta(tx)*Z);
                cB1=CostFunction(H,B1,D);
                if (IL(nn)-cB1)>=beta(tx)*m
                    beta(tx)=2*beta(tx);
                else
                    break;
                end
            end
            
            while beta(tx)>opts.MinStepSize %Step 4
                B2{tx}=retraction(V{tx}+beta(tx)*Z);
                cB2=CostFunction(H,B2,D);
                if (IL(nn)-cB2)<0.5*beta(tx)*m
                    beta(tx)=0.5*beta(tx);
                else
                    break;
                end
            end
        end
        
        V{tx}=retraction(V{tx}+beta(tx)*Z);
    end
    IL(nn+1)=CostFunction(H,V,D);
    
    if mod(nn,opts.Verbose)==0
        fprintf('         Iter.: %3d / %3d,  Interference leakage: %.2e\n',nn,opts.MaxIter,IL(nn+1));
        semilogy(IL(1:nn+1),'b');
        drawnow;
    end
    
    if IL(nn+1)<opts.Tol
        IL=IL(1:nn+1);
        break;
    end
end

%% Calculate ZF decoders
U=cell(1,K);
for rx=1:K
    %Interference covariance matrix
    Q=0;
    for tx=[1:rx-1 rx+1:K]
        Q=Q+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';
    end
    
    %Obtain the subspace that contains the least interference -> Decoders
    [A,temp1,temp2]=svd(Q);
    U{rx}=A(:,end-D(rx,rx)+1:end);  % smallest eigenvectors -> interference free subspace
end

end

function c = CostFunction(H,V,D)
K=length(V);
c=0;
for rx=1:K
    Q=0;
    for tx2=[1:rx-1 rx+1:K]
        Q=Q+H{rx,tx2}*V{tx2}*V{tx2}'*H{rx,tx2}';
    end
    lambdas=svd(Q);
    c=c+sum(lambdas(end-D(rx,rx)+1:end));
end
end


function P=proj_stiefel(X)
[U, S, V] = svd(X, 'econ'); %#ok
P = U*V';
end

function P=proj_grassmann(X)
[P, unused]=qr(X,0);
end

