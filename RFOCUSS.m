function [Xk]=RFOCUSS(A,b,p,iter_max,er)
[~,n]=size(A);
W=eye(n);
X=zeros(n,1);
iter=0;
    while(1)
        Xk=W^2*A'*pinv(A*W^2*A')*b;
        iter=iter+1;
        if norm(Xk-X)/norm(Xk)<er | iter==iter_max
            break
        end
        W=diag(abs(Xk).^(1-p/2));
        X=Xk;
    end
end