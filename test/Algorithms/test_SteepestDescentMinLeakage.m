clc;
clear;
close all;

K = 4;              % Num users
D=diag(2*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system
nT=5*ones(K,1);     % Tx antennas
nR=5*ones(K,1);     % Rx antennas

H = generatechannel(nT,nR,A);

V0=arrayfun(@(a,b)orth(crandn(a,b)),nT',diag(D)','UniformOutput',false);

%%
options.Verbose=100;
options.MaxIter=10000;
options.Tol=1e-8;
options.StartingPoint.V=V0;

options.Manifold='none';%Non-compact Stiefel manifold
[U, V, IL1] = SteepestDescentMinLeakage(H,D,options);

options.Manifold='Stiefel';%Compact Stiefel manifold
[U, V, IL2] = SteepestDescentMinLeakage(H,D,options);

options.Manifold='Grassmann';%Compact Stiefel manifold
[U, V, IL3] = SteepestDescentMinLeakage(H,D,options);

semilogy(IL1,'b');
hold on;
semilogy(IL2,'r');
semilogy(IL3,'k');