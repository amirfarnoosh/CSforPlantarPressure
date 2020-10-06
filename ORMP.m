function [X,er]=ORMP(A,b,k)
[~,n]=size(A);
r=b;
At=[];
index=[];
X=zeros(n,1);
for i=1:k    
    
    J=r'*A/sqrt(sum(r.^2))./sqrt(sum(A.^2));
    %J(abs(J)>1)=0;
    [~,indx]=max(abs(J));
    At=[At A(:,indx)];
    Xt=pinv(At'*At)*At'*b;
    %Xt(i,1)=J(indx)*norm(r)/norm(A(:,indx));
    index=[index;indx];
    r=At*Xt-b;
    
end
X(index)=Xt;
er=r'*r/(b'*b);
end