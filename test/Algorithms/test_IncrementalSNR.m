clc;
clear;
close all;

str='(3x3,1)^5';

[nT nR d K]=system_str2vec(str);

H=generate_channel(nT,nR,K);

for tx=1:K
    for rx=[1:tx-1 tx+1:K]
   H{rx,tx}=1e-1*H{rx,tx};     
    end
end
[U, V] = IncrementalSNR(H);

H{1,1}*V{1}

svd([H{1,2}*V{2} H{1,2}*V{3}])

