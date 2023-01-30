function [Smafra]=runSaFraV1(Y,data,label,dctype)

[L,N]=size(Y);[nr,nc,nb]=size(data);
%----dictionary and regularization parameters
[mfra_alpha,mfra_beta,A]=SaFraPara(dctype);
opts = [];
sal=saliency(Y);
salW=reshape(1./(abs(sal)),N,1)';
%----input
opts.A = A;opts.maxit = 100;opts.sw = salW;opts.c = size(A,2);opts.gt = label;
opts.rho1 = 0.1;opts.rho2 = 0.1;opts.rho3 = 0.1;opts.tol = 1e-4;
opts.frame = 1;opts.Level = 1;  opts.F_it = 1;opts.wLevel = 1/2;opts.x_size = [nr,nc];
opts.X = sunsal(A,Y,'lambda',0,'ADDONE','no','POSITIVITY','yes', 'TOL',1e-4, 'AL_iters',200,'verbose','no');
opts.lambda1 = mfra_alpha;
opts.lambda2 = mfra_beta;
opts.beta1 = 200*opts.lambda2;opts.beta2 = 0.1;
%----run safra
[X,S,Out ] = SaFra(Y,opts);

Smafra=sqrt(sum(S.^2));
Smafra      = (Smafra-min(Smafra(:)))./(max(Smafra(:))-min(Smafra(:)));
end

function [mfra_alpha,mfra_beta,A]=SaFraPara(dctype)

    
if strcmp(dctype,'sandiego')
    load LRMFra_SanDiego_Para
    mfra_alpha=1e-2;mfra_beta=1e-1;
end
    

end

function sal=saliency(Y)
[L,N]=size(Y);
[w Rw] = estNoise(Y);
[kf,Ek,E,delta_p]=hysime(Y,w,Rw,'off');
clear w;
% project in the subspace
p = kf;
Ek = Ek(:,1:p);
Yp = Ek*Ek'*Y;
[Up,D] = svds(Yp*Yp'/N,p);
% project onto the subspace span{E}
sal = sqrt(sum((Y-Up(:,1:1)*Up(:,1:1)'*Y).^2));
sal=(sal-min(sal(:)))/(max(sal(:))-min(sal(:)));
end