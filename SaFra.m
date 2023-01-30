function [X,S,result,it] = SaFra(Y,opts)
%% Initiation
maxit = opts.maxit;
tol = opts.tol;              
A = opts.A;
X = opts.X;
sW = opts.sw;
gt = opts.gt;
rho2 = opts.rho2;
lambda1 = opts.lambda1;    
x_size = opts.x_size;            
frame = opts.frame;           
Level = opts.Level;                                
S = zeros(size(Y)); %Y-A*X;% 
%%
X_p = X;
A_p = A;
Out.Rse = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[D,~] = GenerateFrameletFilter(frame);
m = x_size(1);
n = x_size(2);
e = size(X,1);
X_transition = zeros(m,n*e);
Z = FraDecMultiLevel(X_transition,D,Level);
Theta = Z;
Theta2 = zeros(size(X));
U = zeros(size(X));

result=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it = 1:maxit
    Xup=X;
    Sup=S;
    
    W = ((Y - A * X) + rho2 * S)./(1 + rho2);
    Rw=1./(sqrt(sum(W.^2))+eps);
    S = solve_l1l2(W,sW.*Rw.*lambda1/(1+rho2));

    D = Y - S;
    [Z,Theta,U,Theta2,X] = Framelet_X(A,D,X,opts,Z,Theta,U,Theta2);
     
    %%
    Y1 = A*X;
    Y2 = A_p*X_p;
    Rel_Err = norm(Y1 - Y2,'fro')/norm(Y1,'fro');
    res_S=norm(S-Sup,'inf');
    res_X=norm(X-Xup,'inf');
    
    R1value=sqrt(sum(S.^2));
    [FA1,PD1] = perfcurve(gt,R1value,'1') ;
    auc=-sum((FA1(1:end-1)-FA1(2:end)).*(PD1(2:end)+PD1(1:end-1))/2);

    A_p = A;
    X_p = X;
    if it > 2
    Out.Rse = [Out.Rse;Rel_Err];  
    end
    if Rel_Err < tol  
        break;
    end
    result(it,1)=Rel_Err;result(it,2)=norm(Y-A*X-S,'fro');result(it,3)=res_X;result(it,4)=res_S;result(it,5)=auc;
    
    fprintf('\n[%d] [Rel_Err:%.4f --Obj:%.4f --Xdiff:%.4f --Sdiff:%.4f --AUC:%.4f]\n', it, Rel_Err,norm(Y-A*X-S,'fro'),res_X,res_S,auc);
end

end

function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda(i));
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
end