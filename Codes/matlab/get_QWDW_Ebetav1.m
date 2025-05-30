function [res]=get_QWDW_Ebetav1(alL,t1L,t2L,alR,t1R,t2R,r,k)
beta1=r*exp(1i*k);
[~,~,rL]=get_beta_QW_CS2(alL,t1L,t2L,1);
beta2=rL./beta1;
betaL1=(abs(beta1)>abs(beta2)).*beta1+(abs(beta1)<abs(beta2)).*beta2;
betaL2=(abs(beta1)>abs(beta2)).*beta2+(abs(beta1)<abs(beta2)).*beta1;
if abs(beta1)==abs(beta2)
betaL1=beta1;betaL2=beta2;
end
[MatL,~,~,~]=get_QWDW_mat(t1L,t2L,alL,beta1);
EE=eig(MatL);
[etaR1,etaR2,rR]=get_beta_QW_CS2(alR,t1R,t2R,EE(1));
xy(1)=betaL1.*etaR1;
xy(2)=betaL1.*etaR2;
xy(3)=betaL2.*etaR1;
xy(4)=betaL2.*etaR2;
switch sum(abs(xy)>1)
case 4 
res=abs(xy(4))-1;
case 3
res=abs(xy(4))-1;
case 2
if abs(xy(2))>abs(xy(3))
res=abs(etaR1)-abs(etaR2);
else
res=abs(betaL1)-abs(betaL2);
end
case 1
res=abs(xy(1))-1;
case 0
res=abs(xy(1))-1;
end
%  abs(xy)
end