clc;
clear;
close all;

K = 5;              % Num users
T = 1;
D=diag(1*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system
nT=3*ones(K,1);     % Tx antennas
nR=3*ones(K,1);     % Rx antennas

epsilon=0.1; % Minimum eigenvalue constraint
MaxIter=10; % Maximum number of iterations
IleakTol=1e-10; % Interference leakage tolerance

options.NumExtensions=T;
options.ACS=false;
options.ConstantExtensions=true;
options.A=A;
H = GenerateChannel(nT,nR,options);

RankConstrainedRankMinimization(H,D,epsilon,MaxIter,IleakTol);