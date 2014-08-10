clc;
clear;
close all;

MAXANT=20;
MAXD=5;
MIND=0;
MAXUSERS=5;
MINUSERS=2;


while 1
    M=randi([MINUSERS MAXUSERS]);  %Number of transmitters
    N=randi([MINUSERS MAXUSERS]); %Number of receivers
    D=randi([MIND MAXD],N,M); %Demands matrix
    
    MINTXANT=sum(D)';
    MINRXANT=sum(D,2);
    try
        nT=zeros(M,1);
        for tx=1:M
            nT(tx)=randi([MINTXANT(tx) MAXANT]); %Transmit antennas
        end
        nR=zeros(N,1);
        for rx=1:N
            nR(rx)=randi([MINRXANT(rx) MAXANT]); %Transmit antennas
        end
    catch
        continue;
    end
    
    % K=4;
    % nT=5*ones(K,1);
    % nR=5*ones(K,1);
    % D=diag(2*ones(K,1));
    
    H=GenerateChannel(nT,nR);
    [U V]=RandomBeamforming(H,D);
    B=buildlinearmappingXnetwork_canonical(H,U,V,ones(size(D)),D);
    
    [Ne1, Nv1]=size(B);
    [Ne1 Nv1]
    
    [Ne2, Nv2]=CountEqVarX(nT,nR,D);
    [Ne2 Nv2]
    
    if Ne1~=0 && Ne2~=0 %They may differ is number of equations is zero
    assert(all([Ne1 Nv1]==[Ne2 Nv2]))
    end
end