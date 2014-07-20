function result=weakRICs(n,len,maxPart,options)
%weakRICs: Generates length "len" weak restricted integer compositions (RIC)
% of a number "n" with parts smaller than or equal to "maxPart" (possibly a
% vector)
%
% Reference: inspired by MuPAD's function combinat::compositions
%
% Original source code: 
% sprintf(char(feval(symengine, 'expose', 'combinat::compositions')))

persistent res;
persistent q;
persistent r;

%% Default behavior
opts.BurstLen=Inf; %Default burst length
opts.Append=false; %Maximum number of alternating minimization iterations

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

%%
if not(opts.Append)
   res={}; 
end
% Define reset state
if (isempty(res) && isempty(q) && isempty(r)) || opts.BurstLen==0
res={};  
q=[zeros(1,len-1),n];
r=len;
end

count=0;
while count < opts.BurstLen
    if all(q<=maxPart)
        res{length(res)+1}=q;
        count=count+1;
    end
    if q(len)==0
        if r==1
            q(r)=q(r)+1; %little trick for persistent state
            break;
        else
            q(len)=q(r)-1;
            q(r)=0;
            r=r-1;
        end
    else
        q(len)=q(len)-1;
        r=len-1;
    end
    q(r)=q(r)+1;
end        
result=cell2mat(res');
end