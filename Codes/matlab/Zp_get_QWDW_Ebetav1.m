
% k=pi*2+0.005;
k=0;
t1L= 0.30*pi;
t2L= 0.60*pi;
t1R= 0.60*pi;
t2R= 0.30*pi;
alL=0.3;
alR=0.3;
x0=0.01;
x1=30;
dx=0.01;
for i=1:(x1-x0)/dx
    r(i)=x0+dx*(i-1);
 [res1(i)]=get_QWDW_Ebetav1(alL,t1L,t2L,alR,t1R,t2R,r(i),k);
 [res2(i)]=get_QWDW_Ebetav1(alR,t1R,t2R,alL,t1L,t2L,r(i),k);
end
% figure(11);
% plot(r,res1);hold on;plot(r,zeros(size(r)));
% axis([x0 x1 -1 5])
% figure(21);
% hold on;plot(r,res2);hold on;plot(r,zeros(size(r)));
% axis([x0 x1 -1 5])