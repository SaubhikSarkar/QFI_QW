clear;
t1L= 0.30*pi;
t2L= 0.60*pi;
t1R= 0.60*pi;
t2R= 0.30*pi;
alL=0.3;
alR=0.3;
N=30;
MM=1000;
k=0;
Zp_get_QWDW_Ebetav1;
inibeta=1; %input('Enter the ini1:');
inieta=1; %input('Enter the ini2:');
[k1,wk]=lgwt(MM,-pi,pi);

for i=1:sum(size(k1))-1
    kbeta(i)=k1(i);
keta(i)=k1(i);
if i==1
[r1(i),r2(i),EL(i),checkr(i)]=QWDW_Ebeta2(alL,t1L,t2L,alR,t1R,t2R,kbeta(i),inibeta);
[R1(i),R2(i),ER(i),checkR(i)]=QWDW_Ebeta2(alR,t1R,t2R,alL,t1L,t2L,keta(i),inieta);
else
    
[r1(i),r2(i),EL(i),checkr(i)]=QWDW_Ebeta2(alL,t1L,t2L,alR,t1R,t2R,kbeta(i),abs(r1(i-1)));
[R1(i),R2(i),ER(i),checkR(i)]=QWDW_Ebeta2(alR,t1R,t2R,alL,t1L,t2L,keta(i),abs(R1(i-1)));
end
end
% kk=[-k1;k1(end:-1:1)];
% wkk=[wk;wk(end:-1:1)];
% rr1=[r1';r1(end:-1:1)'];
% rr2=[r2';r2(end:-1:1)'];
% RR1=[R1';R1(end:-1:1)'];
% RR2=[R2';R2(end:-1:1)'];
% WNL1=QW_NHCS_WN_beta1(t1L,t2L,alL,rr1,kk,wkk);
% WNL2=QW_NHCS_WN_beta1(t1L,t2L,alL,rr2,kk,wkk);
% WNR1=QW_NHCS_WN_beta1(t1R,t2R,alR,RR1,kk,wkk);
% WNR2=QW_NHCS_WN_beta1(t1R,t2R,alR,RR2,kk,wkk);


[beta1E,beta2E,eta1E,eta2E,E]=get_QWEj_DW_socon(t1L,t1R,t2L,t2R,alL,alR,N);

figure(1);
% subplot(2,1,1);
hold on;plot(abs(r1).*cos(kbeta),abs(r1).*sin(kbeta),'linewidth',1)
hold on;plot(abs(r1).*cos(-kbeta),abs(r1).*sin(-kbeta),'linewidth',1)
 hold on;plot(r2.*cos(-kbeta),r2.*sin(-kbeta),'linewidth',1)
 hold on;plot(r2.*cos(kbeta),r2.*sin(kbeta),'linewidth',1)
   scatter(real(beta1E),imag(beta1E));hold on;scatter(real(beta2E),imag(beta2E));
 figure(2);
% subplot(2,1,1);
hold on;plot(abs(R1).*cos(keta),abs(R1).*sin(keta),'linewidth',1)
hold on;plot(abs(R1).*cos(-keta),abs(R1).*sin(-keta),'linewidth',1)
 hold on;plot(R2.*cos(-keta),R2.*sin(-keta),'linewidth',1)
  hold on;plot(R2.*cos(keta),R2.*sin(keta),'linewidth',1)
  scatter(real(eta1E),imag(eta1E));hold on;scatter(real(eta2E),imag(eta2E));
% hold on;plot(real(betaL2),imag(betaL2),'linewidth',1)
% figure(2)
% scatter(real(E),imag(E))
% subplot(2,1,2);
% hold on;plot(real(etaR1),imag(etaR1),'linewidth',3);
% hold on;plot(real(etaR2),imag(etaR2),'linewidth',3);