function phat=CN(BL,SL,BR,SR,N,CL,CR)
%Parameters must be column vectors, contrast must be row vector

% 
% ZL = @(BL,SL,N,CL) BL + SL.*CL.^N ;
% ZR = @(BR,SR,N,CR) BR + SR.*CR.^N ;
% PHAT = @(BL,SL,BR,SR,N,CL,CR) cat(3, exp(ZL(BL,SL,N,CL))./(1+exp(ZL(BL,SL,N,CL))+exp(ZR(BR,SR,N,CR))),...
%     exp(ZR(BR,SR,N,CR))./(1+exp(ZL(BL,SL,N,CL))+exp(ZR(BR,SR,N,CR))),...
%     1./(1+exp(ZL(BL,SL,N,CL))+exp(ZR(BR,SR,N,CR))) ) ;

ZL = BL + SL.*CL.^N;
ZR = BR + SR.*CR.^N;
phat = cat(3, exp(ZL)./(1+exp(ZL)+exp(ZR)),...
    exp(ZR)./(1+exp(ZL)+exp(ZR)),...
    1./(1+exp(ZL)+exp(ZR)) ) ;
end