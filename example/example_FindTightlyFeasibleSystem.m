% This scripts finds simple feasible systems where no antenna can be removed
% while maintaining feasibility. This kind of systems are referred to as
% "tightly feasible" systems. Term first introduced in:
% P. Kerret and D. Gesbert, ?Interference alignment with incomplete CSIT sharing,?
% ArXiv preprint available: http://arxiv.org/abs/1211.5380, 2012.

clc;
clear;
close all;

while 1
    %Generate random system of the form ?(M x N,d_i) 
    K=randi([3 10],1);
    nT=randi([2 20],1,1);
    nR=randi([2 20],1,1);
    d=randi([1 min([nT nR])-1],K,1);
    
    [Ne Nv]=CountEqVar(nT*ones(K,1), nR*ones(K,1), d, K);
    
    s=Nv-Ne; %Dimensionality of the solution variety
    if s<min(d) & s>0
       fprintf('s=%d,  (%d x %d, %s)^%d\n',s,nT,nR,mat2str(d'),K);
    end   
end