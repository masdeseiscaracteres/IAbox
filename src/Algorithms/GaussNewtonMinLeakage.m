function [U,V,I,success] = GaussNewtonMinLeakage(H,D,options)
%GAUSSNEWTONMINLEAKAGE  Compute interference alignment solutions for 
%                       interference channels using the Gauss-Newton method
%   [U,V] = GaussNewtonMinLeakage(H,D) computes the solution for a K-user
%   interference channel "H" and demands matrix "D". "H" is K x K cell array 
%   containing all pairwise channel matrices. That is, "H{k,l}" is the MIMO
%   channel matrix from the l-th transmitter to the k-th receiver. Matrix 
%   entries are organized in such a way that "H{k,l}(i,j)" stores the 
%   channel from the j-th transmit antenna in the l-th transmitter to the 
%   i-th receive antenna in the k-th receiver.
%   Input argument "D" is the so-called K x K demands matrix which specifies 
%   the number of streams each transmitter wishes to send to every other 
%   receiver. In particular, "D(k,l)" contains the number of streams the 
%   l-th transmitter wants to send to the k-th receiver. Recall that, for 
%   "D" to describe an interference channel, "D" has to be a diagonal matrix 
%   (i.e. each transmitter only wants to send information to its 
%   corresponding receiver). For the sake of convenience, a length-K vector
%   containing only the diagonal elements of the demands matrix is also accepted.
%       The returned variables "V" and "U" are K x 1 cell arrays containing
%   a set of interference alignment precoders and decoders, respectively; 
%   satisfying the interference alignment conditions. That is,
%
%           U{k}'*H{k,l}*V{l}=0  for all k~=l
%
%   [U,V,I,success] = GaussNewtonMinLeakage(...) also returns the vector
%   "I" containing the evolution of the interference leakage with the
%   iterations of the Gauss-Newton method, and the "success" flag which
%   equals "true" if the predefined tolerance has been reached and "false"
%   otherwise. 
%
%   [U,V,I,success] = GaussNewtonMinLeakage(H,D,options) uses the structure
%   "options" to pass optional arguments to modify the default behavior.
%   Structure "options" can be regarded a list of key-value pairs, i.e. 
%
%          options = struct('Key1',value1,'Key2',value2,...) 
%
%   A list of allowed key-value pairs follows:
%       
%          Key: NwtTol
%          Allowed values: Any positive real number
%          Default value: 1e-20
%          Description: Algorithm stops when the interference leakage is
%                       below NwtTol.
%          
%          Key: MaxNwtIter
%          Allowed values: Any positive integer number
%          Default value: 500
%          Description: Algorithm stops after MaxNwtIter iterations.
%          
%          Key: StartingPoint
%          Allowed values: structure with fields "U" and "V" each storing a
%                          K x 1 cell array with decoders and precoders,
%                          respectively.
%          Default value: a random starting point is generated internally
%          Description: Starting point for the algorithm.
%                       Let say your starting point is U0, V0, then
%                       options.StartingPoint = struct('U',{U0},'V',{V0})
%
%          Key: Verbose
%          Allowed values: Any positive integer number
%          Default value: 0 (no verbosity)
%          Description: Convergence information and plot is provided every
%                       Verbose iterations.
%          
%          Key: Solver
%          Allowed values: 'spqr_solve_tbb','spqr_solve','pinv','auto'
%          Default value: 'auto' (fastest solver is chosen automatically)
%          Description: Solver used to obtain the Gauss-Newton update.
%                       Options 'spqr_solve_tbb' and 'spqr_solve' require
%                       a precompiled MEX function. The SPQR_SOLVE
%                       algorithm exploits the intrinsic sparsity in the 
%                       problem. If it is not available in your system,
%                       the default option chooses Matlab's PINV 
%                       automatically (it will be remarkably slower).
%
%
%   Usage examples:
%
%   Solve the (5x5,2)^4 system:
%       K = 4; %Number of users
%       D = diag(2*ones(K,1));
%       nT = 5; %Number of transmit antennas
%       nR = 5; %Number of receive antennas
%       H = cell(K,K); %Channel's cell array
%       for rx = 1:K
%            for tx = 1:K
%                H{rx,tx}=crandn(nR,nT); %Random channel
%            end
%       end
%       [U,V,I, success] = GaussNewtonMinLeakage(H,D);
%
%   Now, solve the same system with some additional options 
%       options.NwtTol = 1e-10; %Desired interference leakage
%       options.Verbose = 1; %Show progress every n iterations
%       for rx = 1:K
%            U0{rx} = crandn(nR,D(rx,rx)); % Random decoder
%       end
%       for tx = 1:K
%            V0{tx} = crandn(nT,D(tx,tx)); % Random precoder
%       end
%       options.StartingPoint = struct('U',{U0},'V',{V0});
%       [U,V,I, success] = GaussNewtonMinLeakage(H,D,options);
%
% Reference:
% 
% Ó. González, C. Lameiro and I. Santamaría, "A Quadratically Convergent
% Method for Interference Alignment in MIMO Interference Channels,"
% accepted for publication in IEEE Signal Processing Letters, Jul. 2014.

%% Default options
opts.NwtTol=1e-20; %Newton's method tolerance
opts.MaxNwtIter=500; %Maximum number of Newton's method iterations
opts.Verbose=0;
opts.Solver='auto';

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%Choose a solver depending on opts.Solver
switch opts.Solver
    case 'spqr_solve_tbb'
        custom_solve=@(A,b) spqr_solve_tbb(A,b, struct ('solution','min2norm'));
    case 'spqr_solve'
        custom_solve=@(A,b) spqr_solve(A,b, struct ('solution','min2norm'));
    case 'pinv'
        custom_solve=@(A,b) pinv(full(A))*b;
    otherwise %Automatically select the supposedly fastest available solver
        if exist('spqr_solve_tbb','file') %If a MEX function is available for the current architecture
            custom_solve=@(A,b) spqr_solve_tbb(A,b, struct ('solution','min2norm')); %Most efficient SuiteSparseQR minimum 2-norm solution
        elseif exist('spqr_solve','file') %If a MEX function is available for the current architecture
            custom_solve=@(A,b) spqr_solve(A,b,struct('solution','min2norm')); %Most efficient SuiteSparseQR minimum 2-norm solution
        else
            warning('GaussNewton:spqrSolveNotAvailable',...
                ['A compiled version of SPQR_SOLVE could not be found in your system. '...
                'Using PINV instead. ' ...
                'Execution times will be remarkably larger.']);
            custom_solve=@(A,b) pinv(full(A))*b; %Mininum 2-norm slow solution
        end
        
        %Check if this system can run the previously chosen routine
        try
            custom_solve(sparse(1),1);
        catch
            warning('GaussNewton:Undefined',...
                ['Something went wrong.'...
                'Using PINV instead. ' ...
                'Execution times will be remarkably larger.']);
            custom_solve=@(A,b) pinv(full(A))*b;
        end
end

%% Demands and connectivity
% Accept the number of streams in a vector
if min(size(D))==1 %D is vector
    D=diag(D); %convert it into a diagonal matrix
end

% Find network biadjacency matrix, A(i,j)=0 if TX j does not interfere RX i
A=cellfun(@(x) norm(x,'fro')>1e-8,H,'UniformOutput',true);

%% Compute dimensions of precoders and decoders
[N, M]=size(D); %Number of users
nT=cellfun('size',H(1,:),2)'; %Number of transmit antennas
nR=cellfun('size',H(:,1),1); %Number of receive antennas
cT=diag(D); %Number of columns in each precoder
cR=diag(D); %Number of columns in each decoder

%% Generate a random initial point if not passed as an option
if isfield(opts,'StartingPoint')
    U0=opts.StartingPoint.U;
    V0=opts.StartingPoint.V;
else
    U0=arrayfun(@(a,b)orth(crandn(a,b)),nR',cR','UniformOutput',false);
    V0=arrayfun(@(a,b)orth(crandn(a,b)),nT',cT','UniformOutput',false);
end

%% Determine which interfering signals must be cancelled
[lmat kmat]=meshgrid(1:M,1:N);
Dp=arrayfun(@(k,l) sum(D([1:k-1 k+1:end],l)),kmat,lmat);%Dp(i,j) equals the number of interfering streams RX i is receiving from TX j
Arow_part=logical(A.*Dp); %Arow_part(i,j)=true if there are equations involving Hij (i.e. "Hij~=0" AND "interference travelling from TX j to RX i")

%% Compute number of variables and equations
% Nv=sum(D,1)*nT+nR'*sum(D,2); %Number of variables
Nfv=sum(D,1)*(nT-diag(D))+(nR-diag(D))'*sum(D,2); %Number of free variables
Ne=sum(sum(cR*cT'.*Arow_part)); %Number of equations
%% Network connectivity
[rxs txs values]=find(Arow_part); %Set of tuples denoting interfering links (rx,tx)
nlinks=length(values); %Number of active interfering links

%% BEGIN algorithm
[U, V, I, success]=Newton(H,U0,V0,Arow_part,D);
%% END algorithm
return

%% Nested function implementing Newton's corrector
    function [U, V, I, success]=Newton(H,U,V,A,D)
        % Variable initialization
        success = false; %default return state
        qf=@(x) qr(x,0); %orthogonalization function
        I=nan(1,opts.MaxNwtIter); %interference leakage evolution
        
        % Evaluate polynomials (residuals) for the first time
        r=EvalPolynomials(H,U,V);
        I(1)=r'*r;
        
        % Show initial leakage if verbosity enabled
        if opts.Verbose
            fprintf('BEGIN -> Iter.: %3d,  Interference leakage: %.2e\n',0,I(1));
        end
        
        %Check if starting point is a valid solution
        if I(1) < opts.NwtTol
            I=I(1);
            success=true;
            return
        end
        
        %Newton's method iteration start here
        for nn=1:opts.MaxNwtIter
            J=Jacobian(H,U,V,A,D); %Compute Jacobian
            Delta_x=-custom_solve(J,r); %Compute minimum-norm update
            
            % Separate update in parts corresponding to...
            Delta_v=Delta_x(1:nT'*cT);      % precoders...
            Delta_u=Delta_x(nT'*cT+1:end);  % and decoders
            
            % Build increments of precoders and decoders
            Delta_V=cellfun(@reshape,...
                mat2cell(Delta_v,nT.*cT,1),...
                num2cell(nT),num2cell(cT),...
                'UniformOutput',false)';
            Delta_U=cellfun(@(a,b,c)reshape(a,b,c)',...
                mat2cell(Delta_u,nR.*cR,1),...
                num2cell(cR),num2cell(nR),...
                'UniformOutput',false)';
            
            % Update precoders and decoders
            V=cellfun(@plus,V,Delta_V,'UniformOutput',false);
            U=cellfun(@plus,U,Delta_U,'UniformOutput',false);
            
            % Orthogonalization step
            [V,Rv]=cellfun(qf,V,'UniformOutput',false);
            [U,Ru]=cellfun(qf,U,'UniformOutput',false);
            
            % Evaluate residuals
            r=EvalPolynomials(H,U,V);
            
            % Evaluate interference leakage
            I(nn+1)=r'*r;
            
            % Show progress every opts.Verbose iterations
            if mod(nn,opts.Verbose)==0
                fprintf('         Iter.: %3d / %3d,  Interference leakage: %.2e\n',nn,opts.MaxNwtIter,I(nn+1));
                semilogy(I(1:nn+1),'b');
                drawnow;
            end
            
            % Check if requested tolerance is satisfied
            if I(nn+1) < opts.NwtTol
                I=I(1:nn+1);
                success=true;
                break
            end
            
        end
        
        % Show final leakage if verbosity enabled
        if opts.Verbose
            fprintf('END   -> Iter.: %3d,  Interference leakage: %.2e\n',nn,I(end));
            semilogy(0:nn,I,'b');
            grid on;
            xlabel('Iteration number');
            ylabel('Interference leakage');
        end
    end

%% EvalPolynomial nested function
    function r=EvalPolynomials(H,U,V)
        r=zeros(Ne,1);%Homotopy function
        pos=1;
        for kk=1:nlinks %For each interfering link
            rx=rxs(kk); %Current link's receiver
            tx=txs(kk); %Current link's transmitter
            
            temp=U{rx}'*H{rx,tx}*V{tx};
            len=cR(rx)*cT(tx);
            
            r(pos:pos+len-1)=temp(:);
            pos=pos+len;
        end
    end

end

%% Jacobian subfunction
function J = Jacobian(H,U,V,A,D)

[N,M] = size(D);
[nT,cT] = cellfun(@size,V);
[nR,cR] = cellfun(@size,U);

[rxs txs values]=find(A);
nrow_part=length(rxs);

%%Build J_V
ncol_part=M;
J_V_cell=cell(nrow_part,ncol_part);
for link=1:nrow_part
    rx=rxs(link);
    tx=txs(link);
    
    for curr_tx=1:ncol_part
        if curr_tx==tx
            %             DG_V_cell{link,curr_tx} = kron(speye(cT(tx)),U{rx}'*H{rx,tx});
            J_V_cell{link,curr_tx} = IkronA(cT(tx),U{rx}'*H{rx,tx});
        else
            J_V_cell{link,curr_tx} = spalloc(cT(tx)*cR(rx),cT(curr_tx)*nT(curr_tx),0);
        end
    end
end

%%Build J_U
ncol_part=N;
J_U_cell=cell(nrow_part,ncol_part);
for link=1:nrow_part
    rx=rxs(link);
    tx=txs(link);
    
    for curr_rx=1:ncol_part
        if curr_rx==rx
            J_U_cell{link,curr_rx} = kron(V{tx}.'*H{rx,tx}.',speye(cR(rx)));
            %             DG_U_cell{link,curr_rx} = AkronI(V{tx}.'*H{rx,tx}.',cR(rx));
        else
            J_U_cell{link,curr_rx} = spalloc(cT(tx)*cR(rx),nR(curr_rx)*cR(curr_rx),0);
        end
    end
end

%%Concatenate both submatrices
J=cell2mat([J_V_cell J_U_cell]);
end

%% Faster kronecker products when one of the terms is the identity matrix
function X=IkronA(L,A)
[R C] = size(A);
[ia,ja,sa] = find(A);
v=0:L-1;
ix = bsxfun(@plus,R*v,ia(:));
jx = bsxfun(@plus,C*v,ja(:));
X = sparse(ix,jx,sa(:)*ones(1,L),L*R,L*C);
end

function X=AkronI(A,L)
[R C] = size(A);
[ia,ja,sa] = find(A);
v=1:L;
ix = bsxfun(@plus,L*(ia(:)-1).',v(:));
jx = bsxfun(@plus,L*(ja(:)-1).',v(:));
X = sparse(ix,jx,ones(L,1)*A(:).',L*R,L*C);
end
