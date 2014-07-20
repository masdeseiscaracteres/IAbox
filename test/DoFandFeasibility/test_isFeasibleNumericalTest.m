%%
clc;
clear;
close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Test different scenarios                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write a system with the appropriate syntax

str = '(2x2,1)^3';

% Additional examples
% str = '(2x2,1)^3';
% str = '(2x2,1)^4';
% str = '(3x3,2)^2'; %Infeasible according to [3]. It does not satify the outer bound considering the users in pairs. Proper but not feasible
% str = '(4x8,3)^3'; %Proper but not feasible
% str = '(4x8,3)^2(4x8,2)';
% str = '(5x11,4)^3'; %Proper but not feasible
% str = '(6x14,5)^3';
% str = '(6x14,5)(6x14,5)(6x14,4)';
% str = '(5x11,3)^3';
% str = '(8x12,5)^3';
% str = '(8x12,5)(8x12,5)(8x12,4)';
% str = '(4x4,2)^3';
% str = '(5x5,2)^4';
% str = '(2x4,2)^2'; %appeared in [4]
% str = '(3x4,1)^6'; %appeared in [4]
% str = '(4x3,2)^3';
% str = '(2x3,1)^2(3x2,1)^2'; %Proper and Feasible according to [1]
% str = '(2x3,1)^3(3x2,1)^1'; %Feasible
% str = '(4x4,2)(5x3,2)(6x2,2)'; %Not feasible (appeared in [4]) although it satisfies the 2-user channel outer-bound
% str = '(5x11,6)(5x11,3)(5x11,10)'; %Infeasible, violating p2p outer bound
% str = '(6x8,4)(6x8,3)(6x8,3)'; %Interesting, use properness as a bound on the maximum DoF in the whole network
% str = '(5x5,1)(5x5,2)(5x5,2)(5x5,3)';
% str = '(7x13,4)(7x13,5)(7x13,5)';
% str = '(8x12,4)(8x12,5)(8x12,5)';
% str = '(10x18,6)(10x18,7)(10x18,7)';
% str = '(11x17,6)(11x17,7)(11x17,7)';
% str = '(12x16,6)(12x16,7)(12x16,7)';
% str = '(10x19,6)(10x19,7)(10x19,7)'; %Feasible
% str = '(5x11,3)(5x11,4)(5x11,4)';
% str = '(16x8,4)(16x8,6)(16x8,6)'; 
% str = '(16x8,5)(16x8,5)(16x8,6)';
% str = '(4x8,3)(4x8,2)^3'; %Feasible according to [2]
% str = '(6x5,3)^2(6x5,2)';
% str = '(5x6,1)(5x6,2)(5x6,2)(5x6,2)(5x6,2)';
% str = '(2x2,1)(5x5,2)(5x5,2)(8x8,4)'; %[4]'s conjecture, assymetric distribution of antennas achieves higher DoF (vs (5x5,2)^4)
% str = '(3x4,2)(1x3,1)(10x4,2)'; %Proper but not feasible according to [3]. It violates the 2-user cooperative outer-bound (it can be simplified to (4x7,3)(10x4,2))
% str = '(4x7,3)(10x4,2)';
% str = '(2x1,1)(1x2,1)'; %A subset of equations is not proper, hence system is infeasible [3]
% str = '(2x2,1)^3(3x5,1)'; %A subset of equations is not proper, hence system is infeasible [3]
% str = '(3x4,2)(1x3,1)(10x4,2)'; %appeared in [4]. [4] says it is feasible, this test & alternating minimization says it is not
% str = '(4x7,3)(10x4,2)'; %appeared in [4]
% str = '(2x3,1)^3(2x2,1)'; %appeared in [4]
% str = '(2x3,4)(3x3,2)^2';
% str = '(29x11,8)^4';

% Obtain number of TX antennas (nT), number of RX antennas (nR),
% number of streams (d) and users (K) from the string describing the system
[nT nR d K] = system_str2vec(str);

% Specify connectivity matrix
options.A = ones(K); % Assume systems are fully connected
% Choose verbosity level
options.Verbose = 0;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Check feasibility for the given scenario                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Feasible = isFeasibleNumericalTest(nT,nR,diag(d),options);

if Feasible
    fprintf('Scenario %s is FEASIBLE\n',str)
else
    fprintf('Scenario %s is NOT FEASIBLE\n',str)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             REFERENCES                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

[1] C. M. Yetis, "Wireless network capacity," Ph.D., Istanbul Tech.
Univ., Istanbul, 2010 [Online]. Available: http://sites.google.com/site/cenkmyetis/

[2] T. Gou and S. A. Jafar, "Degrees of freedom of the K user M x N MIMO interference channel,"
IEEE Transactions on Information Theory, vol 56, no 12, Dec. 2010

[3] Cenk M. Yetis, Tiangao Gou, Syed A. Jafar and Ahmet H. Kayran, "On Feasibility of Interference Alignment in
MIMO Interference Networks," IEEE Transactions on Signal Processing, vol. 58, no. 9, Sep. 2010

[4] F. Negro, S. P. Shenoy, I. Ghauri, and D. T. M. Slock, "Interference alignment feasibility in constant coefficient MIMO
interference channels," 2010, IEEE Int.l Workshop on Signal Processing Advances on Wireless Communications, (SPAWC)

%}