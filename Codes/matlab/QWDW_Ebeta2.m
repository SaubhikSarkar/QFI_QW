function [r1,r2,E,check1]=QWDW_Ebeta2(alL,t1L,t2L,alR,t1R,t2R,k,ini)

options=optimset('Display','off','MaxFunEvals',2000, ...
    'TolFun',1e-7,'TolX',1e-8);
 [r1,check]=fsolve(@(x)get_QWDW_Ebetav2(alL,t1L,t2L,alR,t1R,t2R,x,k),ini,options);
 check1=0;
if abs(check)>0.001
ini1=ZPpoint_get_QWDW_Ebetav1(alL,t1L,t2L,alR,t1R,t2R,k,ini);
 [r1,check1]=fsolve(@(x)get_QWDW_Ebetav2(alL,t1L,t2L,alR,t1R,t2R,x,k),ini1,options);
if abs(check1)>0.001
    check1=1
    
else
    check1=0;
end
end
[~,~,rL]=get_beta_QW_CS2(alL,t1L,t2L,1);
 r2=rL./r1;

betaL1=r1*exp(1i*k);
[MatL,~,~,~]=get_QWDW_mat(t1L,t2L,alL,betaL1);
EE=eig(MatL);
E=EE(1);
    
% [res]=get_SSHDW_Ebeta1(t1R,t1L,t2,gamL,gamR,k,r1);

end


