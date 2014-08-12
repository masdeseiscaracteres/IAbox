function [nT nR D opts] = GetTestScenario(num)

if isempty(num)
    num=randi(5);
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
        
    case {5,'(5x5,2)^4'}
        % Multi-stream IC: (5x5,2)^4
        K=4;
        D=diag(2*ones(K,1));
        nT=5*ones(K,1);
        nR=5*ones(K,1);
        opts.NumExtensions=1;
        opts.ConstantExtensions=false;
        opts.ACS=false;
        opts.A=ones(K);
      
    otherwise
        error('Undefined scenario');
end