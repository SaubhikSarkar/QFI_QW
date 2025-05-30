function [Mat,Am,Ap,A]=get_QWDW_mat(tt1,tt2,alpp,betta)

s1=[0,1;1,0];
s2=[0,-1i;1i,0];
s3=[1,0;0,-1];
I=[1,0;0,1];

P0=[1,0;0,0];
P1=[0,0;0,1];
M=(expm(s3*alpp));
R1=(expm(-1i*tt1/2*s2));
R2=(expm(-1i*tt2/2*s2));

Fm=R1*P0*R2;
Fp=R1*P1*R2;
Gm=M*R2*P0*R1;
Gp=M*R2*P1*R1;

% A=(Fm*Gp+Fp*Gm);
% Am=(Fm*Gm);
% Ap=(Fp*Gp);

A=[0,0;0,0];
Am=R2*P0*R1*M;
Ap=R2*P1*R1*M;

Mat=A+Ap/betta+Am*betta;

end

