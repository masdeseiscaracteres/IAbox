function str=system_vec2str(nT,nR,d,varargin)

%Convert all arguments into column-vector format
nT=nT(:);
nR=nR(:);
d=d(:);
str=[];

aux1=[nT nR d];

if length(varargin)>=1 && strcmp(varargin{1},'KeepOrder')
    %Do not permute users
    [e,c]=rle(aux1); %Run-length encode the vector idxs
    tab=[e c];
else
    %Permute users (default behaviour)
    [tab,~,loc]=unique(aux1,'rows');
    tab=[tab histc(loc,1:size(tab,1))];
end

for kk=1:size(tab,1)
    str=[str sprintf('(%dx%d,%d)^%d',tab(kk,1),tab(kk,2),tab(kk,3),tab(kk,4))];
end
str=regexprep(str,'(\^1){1}','');

end

function [e,c,idx]=rle(x)
% This function performs Run Length Encoding to the rows of matrix X.
% [E,C,IDX] = RLE(X) returns the element values in E and their number of
% ocurrences in C. IDX is a vector of indexes such that X(IDX,:)=E.

ind=1;
e(ind,:)=x(1,:);
c(ind)=1;
idx(ind)=1;

for i=2:size(x,1)
    if all(x(i-1,:)==x(i,:))
        c(ind)=c(ind)+1;
    else ind=ind+1;
        e(ind,:)=x(i,:);
        c(ind)=1;
        idx(ind)=i;
    end
end
c=c(:);
idx=idx(:);
end









