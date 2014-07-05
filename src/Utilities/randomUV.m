function [U V]=randomUV(nT,nR,d,K)

%% Random set of *-coders
fun=@(n,d) randn(n,d)+1j*randn(n,d);
U = arrayfun(fun,nR,d,'UniformOutput',false);
V = arrayfun(fun,nT,d,'UniformOutput',false);

%% Orthonormal representatives
% fun=@(n,d) qr(randn(n,d)+1j*randn(n,d),0);
% [U, ~] = arrayfun(fun,nR,d,'UniformOutput',false);
% [V, ~] = arrayfun(fun,nT,d,'UniformOutput',false);

end
