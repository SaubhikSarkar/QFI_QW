function ini=ZPpoint_get_QWDW_Ebetav1(alL,t1L,t2L,alR,t1R,t2R,k,ini1)
x0=ini1+1.5;
x1=ini1-0.3;
dx=-0.001;
ini=ini1;
for i=1:(x1-x0)/dx
    r(i)=x0+dx*(i-1);
 [res1(i)]=get_QWDW_Ebetav2(alL,t1L,t2L,alR,t1R,t2R,r(i),k);
if abs(res1(i))<1e-2 ||(i>1 && (res1(i))*res1(i-1)<0) 
   ini=r(i);
   break;
end
end
end