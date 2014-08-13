function [nT nR D opts] = GetTestScenario(num)

if isempty(num)
    num=randi(10);
end

switch num
    % Syntax:
    % case {num,'name'}
    case {1,'(2x2,1)^3'}
        % Single-stream IC: (2x2,1)^3
        K=3;
        D=diag(1*ones(1,K));
        nT=2*ones(K,1);
        nR=2*ones(K,1);
        opts.NumExtensions=1;
        opts.ConstantExtensions=false;
        opts.ACS=false;
        opts.A=ones(K);
        
    case {2,'[(1x1,3)^4,8]'}
        % [(1x1,3)^4,8], proper signaling, time-varying channel extensions, 1.5 DoF
        K=4; % Number of users
        D=diag(3*ones(K,1)); % Streams per user
        nT=1*ones(K,1); % Transmit antennas
        nR=1*ones(K,1); % Receive antennas
        opts.NumExtensions=8; % Number of channel extensions
        opts.ConstantExtensions=false; % If true, constant channel extensions
        opts.ACS=false; % If true use assymetric complex signaling
        opts.A=ones(K); % Connectivity (biadjacency matrix of the graph)
        
    case {3,'[(1x1,1)^4,3]'}
        % [(1x1,1)^4,3], improper signaling, constant channels, 1.33 DoF
        K=4;
        D=diag(2*ones(K,1));
        nT=1*ones(K,1);
        nR=1*ones(K,1);
        opts.NumExtensions=3;
        opts.ConstantExtensions=true;
        opts.ACS=true;
        opts.A=ones(K);
        
    case {4,'[(2x1,2)^3(2x1,3)^3,6]'}
        % [(2x1,2)^3(2x1,3)^3,6], proper signaling, time-varying channels, 2.5 DoF
        K=6;
        D=diag([2 2 2 3 3 3]);
        nT=2*ones(K,1);
        nR=1*ones(K,1);
        opts.NumExtensions=6;
        opts.ConstantExtensions=false;
        opts.ACS=false;
        opts.A=ones(K);
        
    case {5,'[(1x1,2)^3,5]'}
        %[(1x1,2)^3,5], proper signaling, constant channels, 1.2 DoF
        K = 3;
        D=diag(4*ones(K,1));
        nT=1*ones(K,1);
        nR=1*ones(K,1);
        opts.NumExtensions=5;
        opts.ConstantExtensions=true;
        opts.ACS=true;
        opts.A=ones(K);
        
    case {6,'(5x5,2)^4'}
        % Multi-stream IC: (5x5,2)^4
        K=4;
        D=diag(2*ones(K,1));
        nT=5*ones(K,1);
        nR=5*ones(K,1);
        opts.NumExtensions=1;
        opts.ConstantExtensions=false;
        opts.ACS=false;
        opts.A=ones(K);
        
    case {7,'Single-beam IBC'}
        nT=3*ones(3,1);
        nR=4*ones(6,1);
        D=kron(eye(3),[1;1]);
        opts.NumExtensions=1;
        opts.ConstantExtensions=true;
        opts.ACS=false;
        opts.A=ones(size(D));
        
    case {8,'2UserXC_small'}
        % Example from Agustín and Vidal's paper, 5 DoF
        K = 2;
        D = [1 1;2 1]; %Max 5 DoF without channel extensions
        opts.A=ones(K);          % Fully connected system (Xnetwork)
        nT = [4 3]';
        nR = [2 5]';
        opts.NumExtensions=1;
        opts.ConstantExtensions=true;
        opts.ACS=false;
        opts.A=ones(K);
        
    case {9, '2UserXC_large'}
        % Example from Agustín and Vidal's paper, 29 DoF
        K = 2;
        D = [5 8; 5 11];
        nT = [5 8]';     % Tx antennas
        nR = [6 7]';     % Rx antennas
        opts.NumExtensions=3;
        opts.ConstantExtensions=false;
        opts.ACS=false;
        opts.A=ones(K);
        
    case {10,'3UserXN'}
        K=3;
        nT=5*ones(K,1);
        nR=5*ones(K,1);
        D=ones(K);
        opts.NumExtensions=1;
        opts.ConstantExtensions=false;
        opts.ACS=false;
        opts.A=ones(K);
        
    otherwise
        error('Undefined scenario');
end