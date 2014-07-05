function [M N d K]=system_str2vec(str)

%Read system parameters from string
str2=regexprep(str, ')[^\^]', ')^1(');
str_final=regexprep(str2, ')$', ')^1');
aux=sscanf(str_final,'(%dx%d,%d)^%d',[4 Inf]);
Mi=aux(1,:);
Ni=aux(2,:);
di=aux(3,:);
Ki=aux(4,:);

K=sum(Ki);
M=zeros(K,1);
N=zeros(K,1);
d=zeros(K,1);
idx=1;
for ii=1:length(Ki)
    for jj=1:Ki(ii)
      M(idx)=Mi(ii);
       N(idx)=Ni(ii);
       d(idx)=di(ii);
        idx=idx+1;
    end
end