function [B1,B2]=generate_B_circular(M,N)

Dx=sparse(-eye(N));
Dx_1=ones(N-1,1);
Dx_1=sparse(diag(Dx_1,1));
Dx(N,N)=0;
Dx=Dx+Dx_1;
IM=sparse(eye(M));
B1=kron(Dx,IM);%x

Dy=sparse(-eye(M));
Dy_1=ones(M-1,1);
Dy_1=sparse(diag(Dy_1,1));
Dy(M,M)=0;
Dy=Dy+Dy_1;
IN=sparse(eye(N));
B2=kron(IN,Dy);
%B=[B1;B2];
