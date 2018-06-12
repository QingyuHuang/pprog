function [U,B,V] = gkl_bidiag(A)
m=size(A,1);
n=size(A,2);
p=min(m,n);

U=zeros(m,p);
V=zeros(n,p);
B=zeros(p,p);
start=zeros(m,1);
start(1)=1;

beta(1)=norm(start);
U(:,1)=start/beta(1);

alpha(1)=norm(A'*U(:,1));
V(:,1)=(A'*U(:,1))/alpha(1);


for j=1:p-1
    U(:,j+1)=A*V(:,j)-alpha(j)*U(:,j);
    beta(j+1)=norm(U(:,j+1));
    U(:,j+1)=U(:,j+1)/beta(j+1);
    V(:,j+1)=A'*U(:,j+1)-beta(j+1)*V(:,j);
    alpha(j+1)=norm(V(:,j+1));
    V(:,j+1)=V(:,j+1)/alpha(j+1);
end

B=diag(alpha,0)+diag(beta(2:p),-1);
