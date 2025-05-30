
function [beta1,beta2,eta1,eta2,E]=get_QWEj_DW_socon(t1L,t1R,t2L,t2R,alL,alR,N)
s1=[0,1;1,0];
s2=[0,-1i;1i,0];
s3=[1,0;0,-1];
I=[1,0;0,1];
P0=[1,0;0,0];
P1=[0,0;0,1];

syms theta1 theta2 alp

M=[exp(alp),0;0,exp(-alp)];

R1=simplify(expm(-1i*theta1/2*s2));
R2=simplify(expm(-1i*theta2/2*s2));
% Fm=R1*P0*R2;
% Fs=R1*P1*R2;
% Gs=M*R2*P0*R1;
% Gp=M*R2*P1*R1;
Fm=R2*P0*R1;
Fs=R2*P1*R1;
Gs=M;
Gp=M;


t1=[ones(1,N)*t1L,ones(1,N)*t1R];
t2=[ones(1,N)*t2L,ones(1,N)*t2R];
al=[ones(1,N)*alL,ones(1,N)*alR];


for i=1:N*2
    j=i+1-((i+1)>N*2)*N*2;
    l=((i-1)<1)*N*2+i-1;
    A1=eval(subs(Fm,[theta1,theta2,alp],[t1(i),t2(j),al(j)]))*eval(subs(Gp,[theta1,theta2,alp],[t1(i),t2(j),al(j)]));
    A2=eval(subs(Fs,[theta1,theta2,alp],[t1(i),t2(i),al(i)]))*eval(subs(Gs,[theta1,theta2,alp],[t1(i),t2(i),al(i)]));
    Am=eval(subs(Fm,[theta1,theta2,alp],[t1(i),t2(j),al(j)]))*eval(subs(Gs,[theta1,theta2,alp],[t1(j),t2(j),al(j)]));
    Ap=eval(subs(Fs,[theta1,theta2,alp],[t1(j),t2(j),al(j)]))*eval(subs(Gp,[theta1,theta2,alp],[t1(i),t2(j),al(j)]));
%     S((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=A1+A2;
    S((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=0*(A1+A2);
    S((i-1)*2+1:(i-1)*2+2,(i-1)*2+1+2:(i-1)*2+2+2)=Am;
    S((i-1)*2+1+2:(i-1)*2+2+2,(i-1)*2+1:(i-1)*2+2)=Ap;
end
S(1:2,4*N-1:4*N)=S(1+4*N:2+4*N,4*N-1:4*N);
S(4*N-1:4*N,1:2)=S(4*N-1:4*N,1+4*N:2+4*N);
  S=S(1:N*4,1:N*4);
%   [V,D]=eig(S);
%   DD=diag(D);
  syms lam
     DD=eval(vpasolve(det(S-diag(ones(1,4*N))*lam)==0,lam));
  E=1i*log((DD));
[beta1,beta2,rL]=get_beta_QW_CS2(alL,t1L,t2L,DD);
[eta1,eta2,rR]=get_beta_QW_CS2(alR,t1R,t2R,DD);

end







