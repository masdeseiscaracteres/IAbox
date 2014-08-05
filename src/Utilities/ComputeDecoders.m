function U=ComputeDecoders(H,V,D,type,options)

% Some information on how to design decoders:
%
% Vasileios Ntranos and Giuseppe Caire. "A comparison of decoding
% techniques for Interference Alignment." In Communications Control
% and Signal Processing (ISCCSP), 2012 5th International Symposium on,
% pp. 1-5. IEEE, 2012.
%
% H. Sung, S.-H. Park, K.-J. Lee, and I. Lee, "Linear precoder designs for 
% K-user interference channels," IEEE Trans. Wirel. Commun., vol. 9, no. 1,
% pp. 291–301, Jan. 2010.

switch lower(type)
    case 'zf' %Linear zero-forcing of interference
        U=ComputeZFDecoders(H,V,D);
    case 'mmse' %Per-stream SINR-maximizing linear receiver (interference
        %is treated as colored noise)
        U=ComputeMMSEDecoders(H,V,options.W,D);
    otherwise %Random
        
end
end

function U=ComputeZFDecoders(H,V,D)
% ToDo: extend to X networks

% D(k,l)=d_{kl} indicates the k-th receiver is demanding d_{kl} streams from
% the l-th transmitter
[K,L] = size(D);

%Compute the number of interfering streams going from TX l to RX k, Dp(k,l)
[lmat kmat] = meshgrid(1:L,1:K);
Dp = arrayfun(@(k,l) sum(D([1:k-1 k+1:end],l)),kmat,lmat);
Dpl = logical(Dp); %Dpl(k,l) if interference is going from TX l to RX k

[rxs, txs, values] = find(Dpl); %Set of tuples denoting interfering links (rx,tx)

U=cell(K,1);
for rx=1:K %for each receiver
    Q=0; %Receiver covariance matrix
    for tx=find(Dpl(rx,:))
        Q=Q+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';
    end
    %Obtain the subspace that contains the least interference
    [A,unused1,unused2]=svd(Q);
    U{rx}=A(:,end-D(rx,rx)+1:end);  % singular vectors associated to the D(rx,rx) smallest singular values
end
end

function U=ComputeMMSEDecoders(H,V,W,D)
% ToDo: extend to X networks
% W is a Kx1 cell containing the noise covariance matrix at each receiver
% D(k,l)=d_{kl} indicates the k-th receiver is demanding d_{kl} streams from
% the l-th transmitter
[K,L] = size(D);

U=cell(K,1);
for rx=1:K %for each receiver
    Q=W{rx}; %Receiver covariance matrix
    for tx=1:L
        Q=Q+H{rx,tx}*V{tx}*V{tx}'*H{rx,tx}';
    end
    M=H{rx,rx}*V{rx}; %ToDo: extend to X networks
    U{rx}=inv(Q)'*M;
    %ToDo: Constrain to unit-power
end

end