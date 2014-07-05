function R=Rate(H,U,V,W)
% Rate achieved with complex Gaussian signaling and treating interference
% as noise
%
% Based on:
% R. S. Blum, "MIMO capacity with interference," IEEE Journal on Selected 
% Areas in Communications, vol. 21, no. 5, pp. 793–801, Jun. 2003.

% ToDo: some style changes, generalize to X networks

%%----------------------- Description ------------------------------------%
% Obtains the user rates for a MIMO interference channel
% The rate for each user is calculated as
% R(k)=log2(abs(det((eye(N)+Ik)\Qk)), where Ik is the NxN total interference
% matrix and Qk is the transmitted covariance matrix for the k-th user
%
% H: {K,K}(nR,nT) Interference MIMO channel
% U: {K}(nR,d) decoders
% V: {K}(nT,d) precoders
% W: {K}(nR,nR) Noise covariance matrices

% R: Kx1 vector with the user rates
%
% Ó. González, Jan. 2013

V=V(:); %Force V to be a "column cell"
U=U(:); %Force U to be a "column cell"

K=size(H,2); %Number of users
nR=cellfun(@(x) size(x,1),H); %Number of receive antennas
d=cellfun(@(x) size(x,2),U); %Number of streams

R=zeros(K,1);          % store the achievable rates for each user
Qtx_ch=cell(K,K);  % interference covariance matrices (transmitter + channel)  (Size: nR(uu) x nR(uu))
Qrx=cell(K,1); % interference + noise covariance matrices at the receiver (transmitter + channel + receiver) (Size: d(uu) x d(uu))

if isscalar(W)
   Wp=cell(K,1);
   for rx=1:K
       Wp{rx}=W*eye(nR(rx));
   end
else
    Wp=W;
end
for rx=1:K
    for tx=1:K
            % Interference covariance matrix produced by the tx-th Tx into the rx-th Rx%
            Qtx_ch{rx,tx}=H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';  
    end
end

for uu=1:K
    others=1:K;others(uu)=[];   % a vector with all indexes except the uu-th
    Qrx{uu}=sum(cat(3,Qtx_ch{uu,others}),3)+Wp{uu};   
    R(uu)=log2(abs(det(eye(d(uu))+(U{uu}'*Qrx{uu}*U{uu})\(U{uu}'*Qtx_ch{uu,uu}*U{uu}))));
end

%% N.B. Do not confuse with the rate computed assuming MMSE decoding
% for uu=1:K
%     others=1:K;others(uu)=[];   % a vector with all indexes except the uu-th
%     Qrx{uu}=sum(cat(3,Qtx_ch{uu,others}),3)+Wp{uu};   
%     R(uu)=log2(abs(det(eye(nR(uu))+(Qrx{uu})\(Qtx_ch{uu,uu}))));
% end

