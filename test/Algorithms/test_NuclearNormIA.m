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

H = generatechannel(nT,nR,A,struct('NumExtensions',T,'ACS',false,'ConstantExtensions',true));

NuclearNormIA(H,D,epsilon,MaxIter,IleakTol);