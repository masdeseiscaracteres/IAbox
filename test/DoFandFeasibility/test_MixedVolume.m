clc;
clear;
close all;

%% Some example scenarios
% str='(2x2,1)^3';
% str='(4x4,1)^7';
%str='(4x4,2)^3'; %MV=546, 6 solutions
str='(3x5,1)^7';
% str='(5x7,3)^3';

%Some assymmetric scenarios
% str='(2x3,1)(3x2,1)(2x4,1)(2x2,1)'; %Good example
% str='(2x3,1)^2(3x2,1)^2';  %MV=4, feasible
% str='(2x3,1)^2(3x2,1)^2(5x5,1)'; %MV=4, Feasible

%Interesting scenarios
% str='(3x3,2)^2'; %MV=4, but... infeasible!!!
% str='(2x2,1)^3(3x5,1)'; %MV=0, infeasible
% str='(4x3,2)(2x2,1)^2'; %MV=2, it is feasible. 1 solution
% str='(5x11,4)^3'; % MV>0, but infeasible

[nT nR d K]=read_system_str(str);

%% Call the algorithm
% VERBOSE=0 silent mode
% VERBOSE=1 show execution time
% VERBOSE=2 plot each solution graph and pause 0.2 seconds
% VERBOSE=3 plot each solution graph and wait until user presses a key
options.Verbose=1;

[MV S]=IA_MixedVolume_fast(nT,nR,d,options);

%% Show results
if MV>0
    for iSol=1:length(S)
        fprintf('Solution %d of %d\n',iSol,MV)
        fprintf('===================\n');
        disp(S{iSol});
    end
end
