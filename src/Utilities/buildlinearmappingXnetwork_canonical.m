function B = buildlinearmappingXnetwork_canonical(H,U,V,A,D)

M = size(D,2);
N = size(D,1);

getcols=@(x) size(x,1);
nT=cellfun(getcols,V);
nR=cellfun(getcols,U);

[lmat kmat]=meshgrid(1:M,1:N);

sel=@(rx,tx) sum(D(1:rx-1,tx))+1:sum(D(1:rx,tx));
% Dp=arrayfun(@(k,l) sum(D([1:k-1 k+1:end],l)~=0),kmat,lmat);
% Arow_part=and(A,Dp); %Arow_part is 1 if there are equations involving Hkl

Dp=arrayfun(@(k,l) sum(D([1:k-1 k+1:end],l)),kmat,lmat);
Arow_part=(A.*Dp)~=0; %Arow_part is 1 if there are equations involving Hkl

Ap=arrayfun(@(k,l) sum(A([1:k-1 k+1:end],l)~=0),kmat,lmat);
Acol_part=and(D,Ap); %Acol_part(j,p) is 1 if there are equations involving Vjp

[rxs txs ~]=find(Arow_part);
[desired_rxs desired_txs ~]=find(Acol_part);

nrow_part=length(rxs);
ncol_part=length(desired_rxs);

sel=@(rx,tx) sum(D(1:rx-1,tx))+1:sum(D(1:rx,tx));

%Build JV
JVcell=cell(nrow_part,ncol_part);
for kk=1:nrow_part
    rx=rxs(kk);
    tx=txs(kk);
    I=eye(Dp(rx,tx));
    for jj=1:ncol_part
        desired_rx=desired_rxs(jj); %j
        desired_tx=desired_txs(jj); %p
        
        H12=H{rx,tx}(1:sum(D(rx,:)),D(desired_rx,tx)+1:end);
        H22=H{rx,tx}(sum(D(rx,:))+1:end,D(desired_rx,tx)+1:end);
        
        if tx==desired_tx && rx~=desired_rx && (D(desired_rx,tx)~=0) && A(rx,tx)~=0
            Uf=U{rx}(sum(D(rx,:))+1:end,:);
            JVcell{kk,jj} = kron(I(:,1:D(desired_rx,tx)),H12+Uf'*H22); %BUkjl
            I(:,1:D(desired_rx,desired_tx))=[]; %Remove already used columns
        else
            JVcell{kk,jj} = zeros(Dp(rx,tx)*sum(D(rx,:)),D(desired_rx,desired_tx)*(nT(desired_tx)-D(desired_rx,desired_tx)));
        end
    end
end

%Build JU
Acol_part=sum(Arow_part,2)>0;  %Acol_part(p) is 1 if there are equations involving Up

[desired_rxs, ~, ~]=find(Acol_part);

ncol_part=length(desired_rxs);

JUcell=cell(nrow_part,ncol_part);
for kk=1:nrow_part
    rx=rxs(kk);
    tx=txs(kk);
    for jj=1:ncol_part
        desired_rx=desired_rxs(jj);

        if rx==desired_rx && Arow_part(rx,tx)~=0
            
            temp=[];
            for ii=[1:rx-1 rx+1:N]
                H21=H{rx,tx}(sum(D(rx,:))+1:end,1:D(ii,tx));
                H22=H{rx,tx}(sum(D(rx,:))+1:end,D(ii,tx)+1:end);
                V_jl = V{tx}(:,sel(ii,tx));
                Vf=V_jl(size(V_jl,2)+1:end,:);
                temp=[temp; (H21+H22*Vf).'];
            end
            JUcell{kk,jj} = kron(temp,eye(sum(D(rx,:)))); %BVkjl
        else
            JUcell{kk,jj}=zeros(sum(D(rx,:))*Dp(rx,tx),sum(D(desired_rx,:))*(nR(desired_rx)-sum(D(desired_rx,:))));
        end
    end
end

C=[JVcell JUcell];

% Add zeros with compatible dimensions
% for mm=1:size(C,1)
%    for nn=1:size(C,2)
%        if isempty(C{mm,nn})
%            nrows=max(cellfun(@(x)size(x,1),C(mm,:)));
%            ncols=max(cellfun(@(x)size(x,2),C(:,nn)));
%            C{mm,nn}=zeros(nrows,ncols);
%        end
%    end
% end

B=cell2mat(C);
B=sparse(B);

% plot(svd(B));
% drawnow;

% imagesc(abs(B))
% drawnow