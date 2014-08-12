function H = GenerateChannel(nT,nR,options)
% GENERATECHANNEL Generate a multiuser MIMO channel
%
% Input:
%   nT: Lx1 vector with number of transmit antennas at each transmitter
%   nR: Kx1 vector with number of receive antennas at each receiver
%
% Output:
%
%   H: KxL cell array containing random channel matrices
%
% Options:
%
%   NumExtensions: number of channel extensions (default is 1, i.e. no
% channel extension)
%   ConstantExtensions: 1 (default is independent channel extensions)
%   ACS: 1 (default is symmetric complex signaling)

nTxUser=length(nT);
nRxUser=length(nR);

%Default values
opts.A=true(nRxUser,nTxUser);
opts.NumExtensions=1;
opts.ConstantExtensions=false;
opts.ACS=false; %Asymmetric complex signaling

%Overwrite some parameters defined in "opts"
if exist('options','var') && isstruct(options)
    fieldNames=fieldnames(options);
    for iField = 1:length(fieldNames)
        opts.(fieldNames{iField})=options.(fieldNames{iField});
    end
end

A=logical(opts.A);

acsFunc=@(x) [real(x) -imag(x); imag(x) real(x)];

H=cell(nRxUser,nTxUser);
for rx=1:nRxUser
    for tx=1:nTxUser
        %Full-rank channel
        tempH=A(rx,tx)*1/sqrt(2)*(randn(nR(rx),nT(tx))+1j*randn(nR(rx),nT(tx)));
        
        if opts.ACS
            tempH=cell2mat(arrayfun(acsFunc,tempH,'UniformOutput',false));
        end
        
        if opts.ConstantExtensions
            extensionCoeffs=ones(opts.NumExtensions,1);
        elseif opts.ACS
            extensionCoeffs=randn(opts.NumExtensions,1);
        else
            extensionCoeffs=randn(opts.NumExtensions,1)+1j*randn(opts.NumExtensions,1);
        end
        
        H{rx,tx}=kron(diag(extensionCoeffs),tempH);
    end
end

