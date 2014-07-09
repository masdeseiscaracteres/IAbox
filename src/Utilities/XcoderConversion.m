function [W] = XcoderConversion(W,D,conversion_str)

sel=@(rx,tx) sum(D(1:rx-1,tx))+1:sum(D(1:rx,tx));

canonicalize_decoder=@(y) cellfun(@(x)x*inv(x(1:size(x,2),1:size(x,2))),y,'UniformOutput',false);
freedomize_decoder=@(y) cellfun(@(x)x(size(x,2)+1:end,:),y,'UniformOutput',false);
unfreedomize_decoder=@(y) cellfun(@(x) [eye(size(x,2));x],y,'UniformOutput',false);


while ~isempty(conversion_str)
    str=conversion_str(1:2);
    conversion_str(1:2)=[];
    clear Wo;
    switch str(1)
        case 'c' %canonicalize
            switch str(2)
                case 'p'
                    %% Canonicalize precoders
                    [rxs, txs, ~]=find(D);
                    for kk=1:length(txs)
                        rx=rxs(kk); %Current desired receiver
                        tx=txs(kk); %Current transmitter
                        aux=W{tx}(:,sel(rx,tx));
                        Wo{tx}(:,sel(rx,tx))=aux*inv(aux(1:size(aux,2),1:size(aux,2)));
                    end
                    
                case 'd'
                    %% Canonicalize decoders
                    Wo=canonicalize_decoder(W);
            end
        case 'f' %freedomize
            switch str(2)
                case 'p'
                    %% Freedomize precoders
                    [rxs, txs, ~]=find(D);
                    for kk=1:length(txs)
                        rx=rxs(kk); %Current desired receiver
                        tx=txs(kk); %Current transmitter
                        aux=W{tx}(:,sel(rx,tx));
                        Wo{tx}(:,sel(rx,tx))=aux(size(aux,2)+1:end,:);
                    end
                case 'd'
                    %% Freedomize decoders
                    Wo=freedomize_decoder(W);
            end
            
        case 'u' %unfreedomize
            switch str(2)
                case 'p'
                    %% Unfreedomize precoders
                    [rxs, txs, ~]=find(D);
                    for kk=1:length(txs)
                        rx=rxs(kk); %Current desired receiver
                        tx=txs(kk); %Current transmitter
                        aux=W{tx}(:,sel(rx,tx));
                        Wo{tx}(:,sel(rx,tx))=[eye(size(aux,2));aux];
                    end
                case 'd'
                    %% Unfreedomize decoders
                    Wo=unfreedomize_decoder(W);
            end
            
    end
    
    W=Wo;
end
end