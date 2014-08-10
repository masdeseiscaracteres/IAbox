clc;
clear;
close all;

MAXANT=20;
MAXD=6;
MIND=0;
MAXUSERS=7;
MINUSERS=1;
DENSITY=1;

while 1
    M=randi([MINUSERS MAXUSERS]);  %Number of transmitters
    N=randi([MINUSERS MAXUSERS]); %Number of receivers
    D=randi([MIND MAXD],N,M); %Demands matrix
    A=1-logical(sprand(size(D,1),size(D,2),1-DENSITY));
    
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
    
    % Option 1
    [Ne1, Nv1]=CountEqVar(nT,nR,D,struct('A',A));
    [Ne1 Nv1]
    
    % Option 2
    H=GenerateChannel(nT,nR);
    [U V]=RandomBeamforming(H,D);
    B=buildlinearmappingXnetwork_canonical(H,U,V,A,D);
    
    [Ne2, Nv2]=size(B);
    [Ne2 Nv2]
    
    assert(all([Ne1 Nv1]==[Ne2 Nv2]))
end