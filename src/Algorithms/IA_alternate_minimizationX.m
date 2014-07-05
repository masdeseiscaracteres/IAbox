function [sol IL cont_iter]=IA_alternate_minimizationX(H,Deseos)
% Interference Alignment solution via alternate minimization
%
% Inputs:
% H: Channel matrix indexing coefficients as {rx user, tx user}(rx ant, tx ant)
%
% Outputs:
% sol: Struct with fields U and V.
%
% Reference:
% Approaching the capacity of wireless networks through distributed
% interference alignment, Krishna Gomadam, Viveck R. Cadambe, Syed A. Jafar
%
% Pros:
%  - It finds and IA solution
%  - Useful for checking systems feasibility and guessing the number of IA solutions
%
% Cons:
%  - The found IA solution depends on the initialization point
%  - Large number of iterations to converge
%
% Óscar González Fernández
% Advanced Signal Processing Group (GTAS)
% Department of Communications Engineering (DICOM)
% University of Cantabria

sel=@(rx,tx) sum(Deseos(1:rx-1,tx))+1:sum(Deseos(1:rx,tx));

K=length(H);

nT=cellfun('size',{H{1,:}},2);
nR=cellfun('size',{H{:,1}},1);

% if sum(d<=min(nT,nR))<K
% error('Number of antennas at both sides ot the link must be larger than or equal to the number of streams')
% end

Q=cell(1,K);
V=cell(1,K);
U=cell(1,K);

%Alternate minimization algorithm
for tx=1:K
    %Start with arbitrary precoders
    V{tx}=orth(randn(nT(tx),sum(Deseos(:,tx)))+1j*randn(nT(tx),sum(Deseos(:,tx)))); %V'*V=I_d
end

%End conditions
NITER=10e3;
tol=1e-10;

%Algorithm starts here
iter=1;
WLI=nan(1,NITER);
%figure
while true
    for rx=1:K
        %Interference covariance matrix
        Q{rx}=0;
        % for tx=setdiff(1:K,rx)
        for tx=1:K
            for desired_rx=[1:rx-1 rx+1:K]
                Q{rx}=Q{rx}+H{rx,tx}*V{tx}(:,sel(desired_rx,tx))*V{tx}(:,sel(desired_rx,tx))'*H{rx,tx}';
            end
        end
        
        %Obtain the subspace that contains the least interference -> Decoders
        [A,D,B]=svd(Q{rx});
        [eigen,i]=sort(diag(real(D)),'descend');
        A=A(:,i); eigen=eigen(i); % sorted eigenvectors and eigenvalues
        U{rx}=A(:,nR(rx)-sum(Deseos(rx,:))+1:end);  % smallest eigenvectors -> interference free subspace
    end
    
    for tx=1:K
        V{tx}=zeros(nT(tx),sum(Deseos(:,tx)));
        
        for desired_rx=1:K
            %Interference covariance matrix
            Q{desired_rx}=0;
            % for rx=setdiff(1:K,tx)
            for rx=[1:desired_rx-1 desired_rx+1:K]
                Q{desired_rx}=Q{desired_rx}+H{rx,tx}'*U{rx}*U{rx}'*H{rx,tx};
            end
            
            %Obtain the subspace that contains the least interference -> Decoders
            [A,D,B]=svd(Q{desired_rx});
            [eigen,i]=sort(diag(real(D)),'descend');
            A=A(:,i); eigen=eigen(i); % sorted eigenvectors and eigenvalues
            
            V{tx}(:,sel(desired_rx,tx))=A(:,nT(tx)-Deseos(desired_rx,tx)+1:end);  % smallest eigenvectors -> interference free subspace
        end
    end
    
    %Weighted leakage interference
    WLI(iter)=sum(sum(rank_test(Deseos,U,V,H)));
    %   WLI(iter)
    if abs(WLI(iter))<tol
        iter=iter+1;
        break
    end
%     if mod(iter-1,100)==0
%         fprintf('Iter: %d, IL: %g\n',iter,WLI(iter));
% %         semilogy(abs(WLI(1:iter)) );drawnow;
%     end
    iter=iter+1;
end
cont_iter = iter;
WLI=WLI(1:iter-1);
IL=WLI(iter-1);
sol.U=U;
sol.V=V;
%semilogy(abs(WLI) );drawnow;