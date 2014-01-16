function [nT nR d K] = system_str2vec(str)
%SYSTEM_STR2VEC   Obtains a vector representation of an interference
%                 channel specified as a text string.
%   [nT nR d K] = SYSTEM_STR2VEC(str) returns the vectors "nT", "nR" and 
%   "d" containing the number of transmitting antennas, receiving antennas
%   and streams in the system passed in the input string "str". Output
%   argument "K" contains the number of users.
%
%   Usage example:
%      [nT nR d K] = system_str2vec('(3x4,2)(1x3,1)(10x4,2)');

% Copyright (c) 2014, Óscar González Fernández
% All rights reserved.
% Advanced Signal Processing Group (GTAS)
% University of Cantabria (Spain)

%Read system parameters from string
str2=regexprep(str, ')[^\^]', ')^1(');
str_final=regexprep(str2, ')$', ')^1');
aux=sscanf(str_final,'(%dx%d,%d)^%d',[4 Inf]);
nTi=aux(1,:);
nRi=aux(2,:);
di=aux(3,:);
Ki=aux(4,:);

K=sum(Ki);
nT=zeros(K,1);
nR=zeros(K,1);
d=zeros(K,1);
idx=1;
for ii=1:length(Ki)
    for jj=1:Ki(ii)
      nT(idx)=nTi(ii);
       nR(idx)=nRi(ii);
       d(idx)=di(ii);
        idx=idx+1;
    end
end