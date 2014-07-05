function H=inverseIA(U,V,A)

K=length(U);

nR=cellfun('size',U,1);
nT=cellfun('size',V,1);

%% Process connectivity matrix
%Make A optional
if exist('A','var') && ~isempty(A)
    if isequal(size(A),[K K])
        A=(~eye(K)).*A~=0; %Remove diagonal elements, and set non-zeros elements to 1.
    else
        error('Invalid adjacency matrix A');
    end
else
    A=~eye(K);
end

%% Solve inverse IA problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Option 1: solve inverse IA efficiently %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rxs txs values]=find(A); %Set of tuples denoting interfering links (rx,tx)
H=cell(K,K);

% Find orthonormal representatives of *-coders
fun=@(W) qr(W,0);
[Uo, ~] = cellfun(fun,U,'UniformOutput',false);
[Vo, ~] = cellfun(fun,V,'UniformOutput',false);

% Compute W*W' for both U and V
fun=@(W) W*W';
UU_H = cellfun(fun,Uo,'UniformOutput',false);
VV_H = cellfun(fun,Vo,'UniformOutput',false);

% Find a random solution to the inverse IA problem
for kk=1:length(values)
    rx=rxs(kk); %Current link receiver
    tx=txs(kk); %Current link transmitter
    
    %General solution
    X=crand(nR(rx),nT(tx));%Random channel matrix
    H{rx,tx}=X-UU_H{rx}*X*VV_H{tx};   
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Option 2: solve inverse IA using a Kronecker product formulation (less efficient) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rxs txs values]=find(A); %Set of tuples denoting interfering links (rx,tx)
% H=cell(K,K);
% for kk=1:length(values)
%     rx=rxs(kk); %Current link receiver
%     tx=txs(kk); %Current link transmitter
%     Akl=kron(V{tx}.',U{rx}'); %Specify structure of matrices Akl using kronecker products
%     temp_h=null(Akl);
%     temp_h=temp_h*(randn(size(temp_h,2),1)+1j*randn(size(temp_h,2),1)); %Choose a random solution in the nullspace of Akl
%     H{rx,tx}=reshape(temp_h,nR(rx),nT(tx));
% end

end
%% DEBUG: Check interference supression
% for kk=1:length(values)
%     rx=rxs(kk); %Current link receiver
%     tx=txs(kk); %Current link transmitter
%     
%     norm(U{rx}'*H{rx,tx}*V{tx})
% end

%% Complementary information
% Inverse IA with Hermitian channels
% http://www.jstor.org/stable/2100507

% Inverse IA with rank deficient channels
% http://www.jstor.org/stable/25049823

%% Structured inverse IA
% Note: Solving the structured IA problem in some cases requires precoders
% and decoders to have a specific zero pattern
% vec=@(x) x(:);
% 
% M=7;
% N=M;
% d=2;
% 
% H=crandn(N,M);
% mask=1-eye(N,M);
% % mask=1-blkdiag(ones(2,2),ones(2,2),ones(2,2));
% 
% S1=diag(vec(mask));
% S=S1(sum(S1,2)>1e-5,:);
% 
% U=crand(N,d);
% V=crandn(M,d);
% 
% F=orth(U);
% G=orth(V);
% 
% % Suitable for computation
% x=null(S*(eye(N^2)-kron((G*G').',F*F')));
% 
% x=x*crandn(size(x,2),1);%x(:,1);
% X=reshape(x,N,M);
% 
% H_start=X-F*F'*X*G*G';
% 
% % Suitable for theoretical analysis
% S1=diag(vec(~mask));
% S=S1(:,sum(S1,1)>1e-5);
% h=null(kron(V.',U')*S);
% h=h*crandn(size(h,2),1);
% H_start=diag(h);