function [U,V,Cost]=SteepestDescentMinDistance(H,D,options)
% This function is an implementation of the function in
% S. Bazzi, G. Dietl, and W. Utschick, "Interference alignment via
% minimizing projector distances of interfering subspaces," in
% 2012 IEEE 13th International Workshop on Signal Processing Advances in 
% Wireless Communications (SPAWC), 2012, pp. 274–278.
opts.MaxIterAlternating=500;
opts.MaxIterSteepestDescent=100;
opts.minUpdate=1e-4;
opts.Tol=1e-8;
opts.ImprovTol=1e-8;
% opts.MinStepSize=1e-5;

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

Cost=nan(1,opts.MaxIterAlternating);

CostVerbose=cell(K,1);
alpha=zeros(K,1);

P=cell(K,K);
for rx=1:K
    for tx=[1:rx-1 rx+1:K]
        P{rx,tx}=computePkl(H{rx,tx},V{tx});
    end
end

for nn=1:opts.MaxIterAlternating
    gam=ones(K,1); %Initial step size
    for tx=1:K
        %         for tx=randperm(K)
        %Compute local cost function
        alpha(tx)=LocalObjectiveFunction(H,V,P,tx);
        for tt=1:opts.MaxIterSteepestDescent
%             CostVerbose{tx}=[CostVerbose{tx} alpha(tx)];
            
            %Compute derivative with respect to V conjugate
            D_Vconj = Derivative(H,V,P,tx);
            
            %Verify derivative numerically
%             fun=@(X)LocalObjectiveFunction(H,[V(1:tx-1) X V(tx+1:end)],P,tx);
%             [DF_X DF_Xconj]=NumericalDifferentiation(fun,V{tx});
%             assert(norm(D_Vconj-DF_Xconj,'fro')^2<1e-5)
            
            Z=-(eye(nT(tx))-V{tx}*V{tx}')*D_Vconj; %convert from Euclidean gradient to Riemannian gradient
            
            
            %% Armijo's rule
            V_nxt1=V;
            V_nxt2=V;
            m=trace(Z'*Z); %calculate distance on the manifold
            if m<opts.minUpdate
                break;
            end
            while 1
                V_nxt1{tx}=V{tx}+2*gam(tx)*Z;
                nxt_alpha1=LocalObjectiveFunction(H,V_nxt1,P,tx);
                if (alpha(tx)-nxt_alpha1)>=gam(tx)*m
                    gam(tx)=2*gam(tx);
                else
                    break;
                end
            end
            
            while 1%beta(tx)>opts.MinStepSize
                V_nxt2{tx}=V{tx}+gam(tx)*Z;
                nxt_alpha2=LocalObjectiveFunction(H,V_nxt2,P,tx);
                if (alpha(tx)-nxt_alpha2)<0.5*gam(tx)*m
                    gam(tx)=0.5*gam(tx);
                else
                    break;
                end
            end
            
            %% Precoder's update and retraction
            [Q, unused]=qr(V{tx}+gam(tx)*Z,0);
            V{tx}=Q;
            
            alpha(tx)=LocalObjectiveFunction(H,V,P,tx);
        end
        
%         figure(1)
%         semilogy(CostVerbose{tx})
%         hold on;
%         drawnow;
        
        for rx=[1:tx-1 tx+1:K]
            P{rx,tx}=computePkl(H{rx,tx},V{tx});
        end
        
    end
    
    %Calculate global cost function
    Cost(nn)=sum(alpha);
    
    if mod(nn,opts.Verbose)==0
        figure(2)
        fprintf('         Iter.: %3d / %3d,  Cost: %.2e\n',nn,opts.MaxIterAlternating,Cost(nn));
        semilogy(Cost(1:nn),'b');
        drawnow;
    end
  
    if Cost(nn)<opts.Tol || (nn>1 && abs(Cost(nn)-Cost(nn-1))<opts.ImprovTol)
        Cost=Cost(1:nn);
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

function alpha_l = LocalObjectiveFunction(H,V,P,l)
K=size(H,1);
alpha_l=0;
for k=[1:l-1 l+1:K]
    Pkl=H{k,l}*V{l}/(V{l}'*H{k,l}'*H{k,l}*V{l})*V{l}'*H{k,l}';
    for m=setxor(1:K,[l k])
        alpha_l=alpha_l+norm(Pkl-P{k,m},'fro')^2;
    end
end
end

% LocalObjectiveFunction2 is equivalent to LocalObjectiveFunction
% function alpha_l = LocalObjectiveFunction2(H,V,P,l)
% K=size(H,1);
% d=size(V{l},2);
% alpha_l=2*d*(K-1)*(K-2);
% for k=[1:l-1 l+1:K]
%     for m=setxor(1:K,[l k])
%         alpha_l=alpha_l-2*trace(V{l}'*H{k,l}'*P{k,m}*H{k,l}*V{l}*inv(V{l}'*H{k,l}'*H{k,l}*V{l}));
%     end
% end
% alpha_l=real(alpha_l);
% end

function Pkl=computePkl(Hkl,Vl)
Pkl=Hkl*Vl/(Vl'*Hkl'*Hkl*Vl)*Vl'*Hkl';
end

function D_Vconj = Derivative(H,V,P,l)
K=size(H,1);
D_Vconj=0;
for k=[1:l-1 l+1:K]
    %Calculate B
    B=H{k,l}'*H{k,l};
    %Calculate A
    A=0;
    for m=setxor(1:K,[l k])
        A=A+P{k,m};
    end
    A=H{k,l}'*A*H{k,l};
    D_Vconj=D_Vconj+(B*V{l}*inv(V{l}'*B*V{l})*V{l}'*A*V{l}-A*V{l})*inv(V{l}'*B*V{l});
end
D_Vconj=2*D_Vconj;
end
