clc;
clear;
close all;

%% Large system
NSOLS=9;

K=4;
nT=7*ones(K,1);
nR=3*ones(K,1);
D=2*diag(ones(K,1));


H=generatechannel(nT,nR,ones(K));

%% Small system
D2=diag(nR-diag(D));
nT2=(K+1)*diag(D2)-nR;
nR2=nR;
H2=cell(K,K);

for rx=1:K
    for tx=1:K
        H2{rx,tx}=H{rx,tx}(:,1:nT2(tx));
    end
end

%% Calculate solutions for both systems

% Define some options
options.NwtTol=1e-20;
options.MaxNwtIter=1e3;
options.Verbose=0;

%% Run the algorithm for the large system
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
    if length(unique_sols)==NSOLS
        break;
    end
end

sols=unique_sols;
%%
disp('------------------------------------')
%% Run the algorithm for the small system
unique_sols=[];
NSIM=100000;
for nn=1:NSIM
    [U,V,IL,success] = GaussNewtonMinLeakage(H2,D2,options);
    if not(success)
        warning('Solution not found');
        continue;
    end
    %Canonize xcoders
    U=XcoderConversion(U,D2,'cd');
    V=XcoderConversion(V,D2,'cp');
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
    if length(unique_sols)==NSOLS
        break;
    end
end

sols2=unique_sols;


%% Form a contingency table

M=zeros(NSOLS,NSOLS);
for s1=1:NSOLS
    for s2=1:NSOLS
        p=(sols2(s2).U{1})'*sols(s1).U{1};
        M(s1,s2)=norm(p);
    end
end

imagesc(M)


