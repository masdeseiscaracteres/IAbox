function [U V]=UnstructuredApproachLinearBeamforming(H,D)
% Implementation of the USAP-uplink (Unstructured Approach) method described in
%
% Gokul Sridharan and Wei Yu, "Degrees of Freedom of MIMO Cellular 
% Networks: Decomposition and Linear Beamforming Design",
% [Online] Available: http://arxiv.org/abs/1312.2681
%
% Although this method is not proved to work it seems to do so
% for all symmetric proper scenarios where the necessary condition is satisfied
% which makes it an interesting method to take into account
%
% N.B.: It seems to work for systems that have a unique solution
%

[K, K]=size(D); %Number of users
nT=cellfun('size',H(1,:),2)'; %Number of transmit antennas
nR=cellfun('size',H(:,1),1); %Number of receive antennas

%% Assume the network is an IC
d=diag(D);

%% Build matrix M
L=K*d-nR; %number of linear constraints to be satisfied at each receiver
M_cell=cell(K,K);

for rx=1:K
        for tx=1:K
        if tx==rx
            M_cell{rx,tx}=spalloc(L(rx)*nR(rx),d(tx)*nT(tx),0);
        else
            M_cell{rx,tx}=kron(crandn(L(rx),d(tx)),H{rx,tx});
        end
    end
end
M=cell2mat(M_cell);

%% Check necessary condition
% if L'*nR>=nT'*d %Necessary condition is not satisfied (i.e. M is tall)
%     error('The necessary condition L''*nR<nT''*d is not satisfied');
% end
%% Sparse null, it should be equivalent to null(full(M))
[Q R E]=qr(M');
s = abs(diag(R));
tol = norm(M,'fro') * eps(class(M));
r = sum(s > tol);
Z = full(Q(:,r+1:end));

%% Pick a random vector in the nullspace
v=Z*crandn(size(Z,2),1);

%% Reshape v to obtain the sought precoders
V=cellfun(@reshape,...
    mat2cell(v,nT.*d,1),...
    num2cell(nT),num2cell(d),...
    'UniformOutput',false)';

%% Verify precoders are full rank
% dets=cellfun(@(X) det([X crandn(size(X,1),size(X,1)-size(X,2))]),V,'UniformOutput',true);
% if any(abs(dets)<1e-6) %if not throw an error
%     error('A solution could not be found for this scenario.');
% end

%%  Orthonormalization
[V,Rv]=cellfun(@(x) qr(x,0),V,'UniformOutput',false);

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
    U{rx}=A(:,end-d(rx)+1:end);  % smallest eigenvectors -> interference free subspace
end

end


