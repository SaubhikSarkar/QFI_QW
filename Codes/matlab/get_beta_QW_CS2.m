function [betaL1,betaL2,r]=get_beta_QW_CS2(alpha,tt1,tt2,EE)

% beta1= (exp(alpha) + (exp(2*alpha) - 2*EE.^2*exp(2*alpha) + EE.^4*exp(2*alpha) + EE.^2*sin(tt2)^2 + 4*EE.^2*exp(2*alpha)*sin(tt1)^2 + 2*EE.^2*exp(2*alpha)*sin(tt2)^2 + EE.^2*exp(4*alpha)*sin(tt2)^2 + 2*EE.*exp(alpha)*sin(tt1)*sin(tt2) + 2*EE.*exp(3*alpha)*sin(tt1)*sin(tt2) + 2*EE.^3*exp(alpha)*sin(tt1)*sin(tt2) + 2*EE.^3*exp(3*alpha)*sin(tt1)*sin(tt2)).^(1/2) + EE.^2*exp(alpha) + EE.*sin(tt1)*sin(tt2) + EE.*exp(2*alpha)*sin(tt1)*sin(tt2))./(EE.*cos(tt1)*(exp(2*alpha) + cos(tt2) + exp(2*alpha)*cos(tt2) - 1));
% 
% beta2=(exp(alpha) - (exp(2*alpha) - 2*EE.^2*exp(2*alpha) + EE.^4*exp(2*alpha) + EE.^2*sin(tt2)^2 + 4*EE.^2*exp(2*alpha)*sin(tt1)^2 + 2*EE.^2*exp(2*alpha)*sin(tt2)^2 + EE.^2*exp(4*alpha)*sin(tt2)^2 + 2*EE.*exp(alpha)*sin(tt1)*sin(tt2) + 2*EE.*exp(3*alpha)*sin(tt1)*sin(tt2) + 2*EE.^3*exp(alpha)*sin(tt1)*sin(tt2) + 2*EE.^3*exp(3*alpha)*sin(tt1)*sin(tt2)).^(1/2) + EE.^2*exp(alpha) + EE.*sin(tt1)*sin(tt2) + EE.*exp(2*alpha)*sin(tt1)*sin(tt2))./(EE.*cos(tt1)*(exp(2*alpha) + cos(tt2) + exp(2*alpha)*cos(tt2) - 1));
% 
% r=(cos(tt2) - exp(2*alpha) + exp(2*alpha)*cos(tt2) + 1)/(exp(2*alpha) + cos(tt2) + exp(2*alpha)*cos(tt2) - 1);

s1=[0,1;1,0]; s2=[0,-1i;1i,0]; s3=[1,0;0,-1]; I=[1,0;0,1];
P0=[1,0;0,0]; P1=[0,0;0,1];
M=(expm(s3*alpha));
R1=(expm(-1i*tt1/2*s2)); R2=(expm(-1i*tt2/2*s2));

A=[0,0;0,0];
Am=R2*P0*R1*M;
Ap=R2*P1*R1*M;

a1 = Am(1, 1); a2 = Am(1, 2);
a3 = Am(2, 1); a4 = Am(2, 2);
b1 = Ap(1, 1); b2 = Ap(1, 2);
b3 = Ap(2, 1); b4 = Ap(2, 2);
c1 = 0; c2 = 0; c3 = 0; c4 = 0;

a = a1.*c4 + a4.*c1 - a2.*c3 - a3.*c2 - EE.*(a1 + a4 - a2 - a3);
b = a1.*b4 + a4.*b1 - a2.*b3 - a3.*b2 - EE.*(c1 + c4 - c2 - c3);
c = b1.*c4 + b4.*c1 - b2.*b3 - b3.*c2 - EE.*(b1 + b4 - b2 - b3);   
    
beta1 = (-b + (b.^2 - 4*a.*c).^(1/2)) ./ (2*a);
beta2 = (-b - (b.^2 - 4*a.*c).^(1/2)) ./ (2*a);
r = (abs((exp(2*alpha)*sin(tt1/2)*sin(tt2/2) - cos(tt1/2)*cos(tt2/2))./ ...
         (sin(tt1/2)*sin(tt2/2) - exp(2*alpha)*cos(tt1/2)*cos(tt2/2))))^(1/2);

betaL1=(abs(beta1)>abs(beta2)).*beta1+(abs(beta1)<abs(beta2)).*beta2+(abs(beta1)==abs(beta2)).*beta1;

betaL2=(abs(beta1)>abs(beta2)).*beta2+(abs(beta1)<abs(beta2)).*beta1+(abs(beta1)==abs(beta2)).*beta2;

end
