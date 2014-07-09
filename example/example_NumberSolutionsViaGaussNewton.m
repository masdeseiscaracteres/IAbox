%This script finds distinct solutions for a given interference channel using the
%GaussNewton method
clc;
clear;
close all;

%% Some example scenarios
%% 
K = 6;
nT=3*ones(K,1);
nR=11*ones(K,1);
D=diag(2*ones(K,1));

K = 6;
nT=5*ones(K,1);
nR=23*ones(K,1);
D=diag(4*ones(K,1));

K = 4;
nT=6*ones(K,1);
nR=9*ones(K,1);
D=diag(3*ones(K,1));

K = 4;
D=diag(1*ones(K,1));
A=ones(K);
nT=[2 2 3 3];
nR=[3 3 2 2];

K = 3;
nT=7*ones(K,1);
nR=8*ones(K,1);
D=diag([3 3 5]);

K = 4;
nT=3*ones(K,1);
nR=3*ones(K,1);
D=diag(1*ones(K,1));
%% (4x3,2)(2x2,1)^2
K=3;
nT=[4 2 2];
nR=[3 2 2];
D=diag([2 1 1]);
%% (2x3,1)(3x2,1)(2x4,1)(2x2,1)
% K=4;
% nT=[2 3 2 2];
% nR=[3 2 4 2];
% D=diag(ones(K,1));
%% (5x5,2)^4
% K = 4;               % Number of users
% nT=5*ones(K,1);      % Tx antennas
% nR=5*ones(K,1);      % Rx antennas
% D=diag(2*ones(K,1)); % Demands matrix
%% (3x3,1)^5
% K = 5;
% nT=3*ones(K,1);
% nR=3*ones(K,1);
% D=diag(1*ones(K,1));
%% (10x10,1)^19
% K = 19;
% nT=10*ones(K,1);
% nR=10*ones(K,1);
% D=diag(1*ones(K,1));
%% (12x12,4)^5
% K = 5;
% nT=12*ones(K,1);
% nR=12*ones(K,1);
% D=diag(4*ones(K,1));
%% (29x11,8)^4
% K = 4;
% nT=29*ones(K,1);
% nR=11*ones(K,1);
% D=diag(8*ones(K,1));
%% (2x3,1)^2(3x2,1)^2
% K=4;
% nT=[3 3 2 2];
% nR=[2 2 3 3];
% D=diag([1 1 1 1]);
%% (4x4,2)^3
% K=3;
% nT=4*ones(K,1);
% nR=4*ones(K,1);
% D=diag(nT/2);
%% Prepare input arguments
% Generate a random interference channel
H = cell(K,K); %Channel's cell array
for rx = 1:K
    for tx = 1:K
        H{rx,tx}=crandn(nR(rx),nT(tx)); %Random channel matrix
    end
end
% Define some options
options.NwtTol=1e-20;
options.MaxNwtIter=1e3;
options.Verbose=0;

%% Now run the algorithm
unique_sols=[];
NSIM=100000;
for nn=1:NSIM
    [U,V,IL,success] = GaussNewtonMinLeakage(H,D,options);
    if not(success)
        warning('Solution not found');
        continue;
    end
    %Canonize xcoders
    U=XcoderConversion(U,D,'cd');
    V=XcoderConversion(V,D,'cp');
    nsols=length(unique_sols);
    if nsols==0
        unique_sols.U=U;
        unique_sols.V=V;
    else
        distinct=zeros(1,nsols);
        for kk=1:nsols
            dif=sum([cellfun(@(x,y) norm(x-y,'fro')^2,unique_sols(kk).U,U) cellfun(@(x,y) norm(x-y,'fro')^2,unique_sols(kk).V,V)]);
            if dif>1e-3
                distinct(kk)=1;
            else
                break;
            end
        end
        if all(distinct)
            %Append a new solution
            unique_sols(nsols+1).U=U;
            unique_sols(nsols+1).V=V;
        end
    end
     fprintf('Found %d distinct solutions out of %d\n',length(unique_sols),nn);
     
end