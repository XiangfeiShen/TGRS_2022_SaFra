function [Smafra]=runSaFraV1(Y,data,label,dctype)

[L,N]=size(Y);[nr,nc,nb]=size(data);
[mfra_alpha,mfra_beta,A]=SaFraPara(dctype);
opts = [];
saltype='proj2';
sal=shenConsSal(data,100,saltype,'yes','yes');
salW=reshape(1./(abs(sal)),N,1)';
opts.A = A;opts.maxit = 150;opts.sw = salW;opts.c = size(A,2);opts.gt = label;
opts.rho1 = 0.1;opts.rho2 = 0.1;opts.rho3 = 0.1;opts.tol = 1e-4;
opts.frame = 1;opts.Level = 1;  opts.F_it = 1;opts.wLevel = 1/2;opts.x_size = [nr,nc];
opts.X = sunsal(A,Y,'lambda',0,'ADDONE','no','POSITIVITY','yes', 'TOL',1e-4, 'AL_iters',200,'verbose','no');
opts.lambda1 = mfra_alpha;
opts.lambda2 = mfra_beta;
opts.beta1 = 200*opts.lambda2;opts.beta2 = 0.1;

[X,S,Out ] = LRNMF_HSI_AD_RW(Y,opts);

Smafra=sqrt(sum(S.^2));
Smafra      = (Smafra-min(Smafra(:)))./(max(Smafra(:))-min(Smafra(:)));
end