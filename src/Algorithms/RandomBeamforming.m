function [U V]=RandomBeamforming(H,D,options)

% [~, K]=size(H);
nT=cellfun('size',H(1,:),2)';
nR=cellfun('size',H(:,1),1);

opts.ACS=false; %Default is symmetric complex signaling

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

if opts.ACS
    custom_randn=@(r,c)randn(r,c);
else
    custom_randn=@(r,c)randn(r,c)+1j*randn(r,c);
end

Num_Tx = length(D(1,:));
Num_Rx = length(D(:,1));

sel=@(rx,tx) sum(D(1:rx-1,tx))+1:sum(D(1:rx,tx));

% Generate random precoders and decoders

U = cell(Num_Rx,1);
V = cell(Num_Tx,1);

for tx = 1:Num_Tx
    for rx = 1:Num_Rx
        [auxV,~] = qr(custom_randn(nT(tx),sum(D(rx,tx))));
        V{tx}(:,sel(rx,tx)) = auxV(:,1:sum(D(rx,tx)));
    end
end

for rx = 1:Num_Rx
    [auxU,~] = qr(custom_randn(nR(rx),sum(D(rx,:))));
    U{rx} = auxU(:,1:sum(D(rx,:)));
end

end
