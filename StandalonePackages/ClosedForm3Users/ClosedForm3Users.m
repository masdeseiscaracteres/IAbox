function sols=ClosedForm3users(H,idx_sol)
% Interference Alignment closed form solution for the case
%   (2dx2d,d)^3
%       2d transmit antennas
%       2d receive antennas
%       d streams per user (d degrees of freedom, d DoFs)
%       3 users
%
% Inputs:
% H: Channel matrix indexing coefficients as {rx user, tx user}(rx ant, tx ant)
% idx_sol: Indexes of desired solutions (larger eigvalue first). If empty all
% solutions are returned.
%
% Outputs:
% sols: Struct array of size 1xlength(idx_sol) with fields U and V.
%
% Reference:
% V. R. Cadambe and S. A. Jafar, "Interference alignment and degrees of 
% freedom of the K-User interference channel,” IEEE Trans. Inf. Theory, 
% vol. 54, no. 8, pp. 3425–3441, Aug. 2008.

%TODO
% Normalize to distinguish different solutions at a glance

K=3; %Users
[N M]=size(H{1,1}); %Transmitter & Receiver antennas
d=M/2; %Degrees of freedom/user

%% Error checking
if ~isequal(size(H),[3 3])
    error('H is not a 3-user channel')
end
if (mod(M,2)==1) & (mod(N,2)==1)
    error('Number of transmitting and receiving antennas must be even');
end
if M~=N
    error('Number of transmitting and receiving antennas are unequal');
end


%% Closed form solution
mv1=inv(H{3,1})*H{3,2}*inv(H{1,2})*H{1,3}*inv(H{2,3})*H{2,1};
[v1 dv1]=eig(mv1);
[aux idx_dv1]=sort(abs(diag(dv1)),'descend');

cmbs=combnk(1:2*d,d);
if nargin==2
    if sum(idx_sol>size(cmbs,1))>0
        error(sprintf('There are %d different solutions, check idx_sol',size(cmbs,1)));
    end
    cmbs=cmbs(idx_sol,:);
end

for kk=1:size(cmbs,1)
    V{1}=v1(:,idx_dv1(cmbs(kk,:)));
    aux=inv(H{3,2})*H{3,1}*V{1};
    V{2}=aux/norm(aux);
    
    aux=inv(H{2,3})*H{2,1}*V{1};
    V{3}=aux/norm(aux);
    
    [u1 du1]=eig([H{1,2}*V{2} H{1,3}*V{3}]');
    [aux idx]=sort(abs(diag(du1)),'ascend');
    U{1}=u1(:,idx(1:d));
    
%     U{1}'*[H{1,2}*V{2} H{1,3}*V{3}] %0 if aligned
    
    
    [u2 du2]=eig([H{2,1}*V{1} H{2,3}*V{3}]');
    [aux idx]=sort(abs(diag(du2)),'ascend');
    U{2}=u2(:,idx(1:d));
    
%      U{2}'*[H{2,1}*V{1} H{2,3}*V{3}] %0 if aligned
    
    
    [u3 du3]=eig([H{3,2}*V{2} H{3,1}*V{1}]');
    [aux idx]=sort(abs(diag(du3)),'ascend');
    U{3}=u3(:,idx(1:d));
    
    % Verify perfect alignation
    % WLI=calc_WLI(U,H,V,K);

    % Normalize and save output
    % TODO: normalize as [Jafar2009, Feasibility Conditions for Interference Alignment]
    sols(kk).U=cellfun(@(v)v/norm(v),U,'UniformOutput',false);
    sols(kk).V=cellfun(@(v)v/norm(v),V,'UniformOutput',false);
end