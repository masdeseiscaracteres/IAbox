function [MV S]=MixedVolume(nT,nR,D,options)
%
%   MV=MixedVolume(nT,nR,D,options)
%
% Calculates mixed volume of the polynomial equation system
% that describes a general interference alignment (IA) scenario.
% The mixed volume matches the exact number of solutions for all 1-stream
% cases, otherwise it is a bound of the number of solutions.
%
% options.Adj is the biadjacency matrix of the bipartite graph that describes
% network connectivity. The adjacency matrix can be obtained as
% A=[0 B; B' 0] where B is the biadjacency matrix. If it is not passed as
% an argument, the network is assumed to be fully connected.
%
% Example: Scenario (2x2,1)^3
% -> Inputs:
% nT = [2 2 2]; %number of TX antennas
% nR = [2 2 2]; %number of RX antennas
% d  = [1 1 1]; %number of DoF
% [MV S]=MixedVolume(nT,nR,d)
%
% -> Outputs:
% MV =              %mixed volume
%      2
% S =               %solutions set
%     [3x3 double]    [3x3 double]
% S{:}
% ans =
%    NaN     1     0
%      0   NaN     1
%      1     0   NaN
% ans =
%    NaN     0     1
%      1   NaN     0
%      0     1   NaN
%
% References:
% Ó. González, C. Beltrán, and I. Santamaría, "On the Number of 
% Interference Alignment Solutions for the K-User MIMO Channel 
% with Constant Coefficients,” ArXiv preprint available: http://arxiv.org/abs/1301.6196v2, Jan. 2013.

%% Default options
K=length(nT);
opts.Adj=~eye(K); %Adjacency matrix of a fully connected scenario
opts.Verbose=0; % VERBOSE=0 silent mode
                % VERBOSE=1 show execution time
                % VERBOSE=2 plot each solution graph and pause 0.2 seconds
                % VERBOSE=3 plot the solution graph and wait until user press a key
                
%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%%
nT=nT(:);
nR=nR(:);
if isvector(D)
    d=D(:);
else
    d=diag(D);
end
%%
int_graph=opts.Adj;

%Variables initialization
start_tx=1;
start_rx=1;
MV=0;
S=cell(0);
sol_set_dim=NaN;

if K~=size(int_graph,1) || K~=size(int_graph,2)
    error('Number of users and interference graph matrix dimensions must agree');
end

[Ne Nv]=countEqVar(nT, nR, d, K,int_graph);

if Nv<Ne
    %     error('System is improper');
    return;
end

%Dimensions of table
N=sum(d);
int_graph_aux=zeros(N);

%Expand the interference graph to deal with multiple streams
indexes = cumsum([1; d]);
[r c]=find(int_graph);
for kk = 1:length(r)
    int_graph_aux(indexes(r(kk)):indexes(r(kk)+1)-1,indexes(c(kk)):indexes(c(kk)+1)-1) = 1;
    
end
int_graph=int_graph_aux;

%Number of available dimensions pero TX (row) and RX (col)
avail_dimTX_s=nT-d; %available dim per TX and stream
avail_dimRX_s=nR-d; %available dim per RX and stream

avail_dimTX=zeros(1,N);
avail_dimRX=zeros(1,N);

ii=1;
for kk=1:K
    for kkk=1:d(kk)
        avail_dimTX(ii)=avail_dimTX_s(kk); %available dim per TX
        avail_dimRX(ii)=avail_dimRX_s(kk); %available dim per RX
        ii=ii+1;
    end
end

%Plot graph
switch opts.Verbose
    case 1
        %Start timer
        tstart=tic;
        tic;
    case {2,3}
        plot_bipartite_graph(int_graph);
end

%Call backtracking algorithm
[MV S sol_set_dim]=backtracking(MV,S,sol_set_dim,int_graph,avail_dimTX,avail_dimRX,start_tx,start_rx,opts);
if sol_set_dim~=0
    fprintf('Solution set dimensionality larger than or equal to %d\n',sol_set_dim);
end

switch opts.Verbose
    case 1
    %Stop timer
    toc(tstart);
end

end

%%
function [MV S sol_set_dim]=backtracking(MV,S,sol_set_dim,int_graph,avail_dimTX,avail_dimRX,start_tx,start_rx,opts)

if sum(sum(int_graph)<=avail_dimRX)==size(int_graph,1)
    
    %Available dimensions at the RX
    rx_free_vars=sum(avail_dimRX-sum(int_graph));
    %Available dimensions at the RX
    tx_free_vars=sum(avail_dimTX);
    
    total_free_vars=rx_free_vars+tx_free_vars;
    
    if total_free_vars>0
        sol_set_dim=total_free_vars; %solution set dimensionality;
        MV=Inf;
        S=cell(0);
        return;
    end
    
    MV=MV+1;
    S{MV}=int_graph;
    sol_set_dim=0;
    
    %Plot graph
    switch opts.Verbose
        case 1
            %Show elapsed time every 1000 solutions
            if mod(MV,1000)==0
                fprintf('Solutions [%d,%d]\n',MV-1000+1,MV)
                toc;
                tic;
            end
        case {2,3}
            fprintf('Solution %d\n',MV)
            fprintf('===============\n');
            disp(int_graph);
            opts.Type='sol';
            h=plot_bipartite_graph(int_graph,opts);
            drawnow;
            delete(h);
    end
    
    
    
    
else
    
    [cand_tx rx]=construct_candidates(int_graph,avail_dimTX,avail_dimRX,start_tx,start_rx);
    
    for kk=1:length(cand_tx)
        
        %Recover variables from stack
        next_int_graph=int_graph;
        next_avail_dimTX=avail_dimTX;
        next_avail_dimRX=avail_dimRX;
        
        %Spend a dimension of the TX cand_tx(kk)
        next_avail_dimTX(cand_tx(kk))=next_avail_dimTX(cand_tx(kk))-1; %spend a dimension
        %to null a link from tx to rx
        next_int_graph(cand_tx(kk),rx)=next_int_graph(cand_tx(kk),rx)-1; %remove a link
%         next_int_graph
        
%         d=dbstack;
%         sum(cellfun(@(x) strcmp(x,'backtracking'),{d(:).name}))
        
        %Define next tx and
        if sum(next_int_graph(:,rx))==next_avail_dimRX(rx)
            next_start_rx=rx+1;
            next_start_tx=1;
        else
            next_start_rx=rx;
            next_start_tx=cand_tx(kk);
        end
        
        %Recursive call
        [MV S sol_set_dim]=backtracking(MV,S,sol_set_dim,next_int_graph,next_avail_dimTX,next_avail_dimRX,next_start_tx,next_start_rx,opts);
    end
    
end

end

%%
function [cand_tx cand_rx]=construct_candidates(int_graph,avail_dimTX,avail_dimRX,start_tx,start_rx)

N=size(int_graph,1);
%Initialize the array of candidates
cand_tx=[];
cand_rx=[];

% If a RX has an affordable number of interferers
while(sum(int_graph(:,start_rx))<=avail_dimRX(start_rx))
    start_rx=start_rx+1;
end
if start_rx>N
    return;
end

cand_rx=start_rx;
kk=1;
for tx = start_tx:N %TODO:limit the number of candidates to the reasonable minimum
    if avail_dimTX(tx)>0 & int_graph(tx,cand_rx)>0%if this TX can null a not-nulled link
        cand_tx(kk)=tx; %add TX to the list of candidates
        kk=kk+1;
    end
end
end

%%
function h=plot_bipartite_graph(biadj_matrix,options)

opts.Verbose=0;
opts.Type='';

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

h=[];
K=size(biadj_matrix,1);

prop=3/4;
tx=-prop*K/2;
rx=prop*K/2;
user1=(K-1)/2;
userK=-(K-1)/2;


for t=1:size(biadj_matrix,1)
    for r=1:size(biadj_matrix,2)
        
        if biadj_matrix(t,r)>0
            
            if strcmp('sol',opts.Type)
                h=[h plot([tx rx],user1+1-[t r],'Color',[0 0 0],'LineWidth',2)];
            else
                h=[h plot([tx rx],user1+1-[t r],':','Color',[0.8 0.8 0.8])];
            end
            hold on;
        end
    end
end

for t=1:size(biadj_matrix,1)
    for r=1:size(biadj_matrix,2)
        
        h=[h plot([tx rx],user1+1-[t r],'bo','MarkerSize',7,'MarkerFaceColor','b','LineWidth',2)];
        h=[h text(1.2*tx,user1+1-t,sprintf('T_{%d}',t), ...
            'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle')];
        h=[h text(1.2*rx,user1+1-r,sprintf('R_{%d}',r), ...
            'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle')];
    end
end

axis([1.1*[tx rx] 1.2*[userK user1]]);
axis equal
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
switch opts.Verbose
    case 2
        pause(0.2);
    case 3
        pause;
end
end
