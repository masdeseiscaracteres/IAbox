function I = InterferenceLeakage(H,U,V,D)
% I = InterferenceLeakage(H,U,V,D)
% H: cell array containing MIMO channels
% U: cell array containing decoders
% V: cell array containing precoders
% D: demands matrix

Num_Tx = length(D(1,:));
Num_Rx = length(D(:,1));

I = zeros(Num_Rx,Num_Tx);

sel=@(rx,tx) sum(D(1:rx-1,tx))+1:sum(D(1:rx,tx));

% Evaluate interference leakage
for rx = 1:Num_Rx
    for tx = 1:Num_Tx
        V_int = V{tx};
        V_int(:,sel(rx,tx)) = [];
        I(rx,tx) = norm(U{rx}'*H{rx,tx}*V_int, 'fro')^2;
    end
end