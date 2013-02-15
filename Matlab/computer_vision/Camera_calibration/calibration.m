%%
% CS 6320 : 3D Computer Vision
%
% Author : Arthur COSTE
% Date : January 2013
%
% Content : Camera calibration
%
%%

clear all
close all

I=imread('P1020660.JPG');
I2=double(I(:,:,1));

imagesc(I2)

%[x,y] = select_points(I2)

%small image
%x= [261.9035;335.5877;398.7456;461.9035;514.5351;602.2544;651.3772;696.9912;
  %760.1491;816.2895;889.9737;258.3947;335.5877;398.7456;458.3947;511.0263;
  %605.7632;647.8684;700.5000;753.1316;823.3070;893.4825;893.4825]

%y= [743.9211;730.7632;717.6053;701.8158;688.6579;691.2895;707.0789;722.8684;
  %741.2895;759.7105;783.3947;667.6053;657.0789;646.5526;638.6579;628.1316;
  %630.7632;638.6579;651.8158;662.3421;678.1316;696.5526;693.9211]
  
  %original image
%x= 1.0e+03 .*[0.7845;1.0192;1.2218;1.4032;1.5738;1.8405;1.9792;2.1178;2.2992;
    %2.4805;2.7152;0.7738;0.9978;1.2112;1.3925;1.5525;1.8298;1.9685;2.1178;
    %2.2992;2.4912;2.7152;3.3018]

%y = 1.0e+03.*[2.2605;2.2205;2.1725;2.1325;2.0925;2.1005;2.1405;2.1885;2.2525;
    %2.3005;2.3805;2.0285;2.0045;1.9645;1.9325;1.9005;1.9005;1.9405;1.9725;
    %2.0205;2.0605;2.1165;2.0205]
    
x = 1.0e+03.*[0.7845;1.0192;1.2112;1.4032;1.5632;1.8192;1.9792;2.1285;2.2992;
    2.4805;2.7258;0.7845;1.0085;1.2005;1.4032;1.5525;1.8298;1.9685;2.1285;
    2.2885;2.4912;2.7045;0.7845;0.7845;0.7738;0.7738;0.7632;0.7845;2.7258;
    2.7152;2.7152;2.7258;2.7152;2.7152;3.2058]

y = 1.0e+03.*[2.2605;2.2205;2.1725;2.1325;2.1005;2.1005;2.1485;2.1885;2.2525;
    2.3005;2.3885;2.0285;1.9885;1.9565;1.9325;1.9005;1.9005;1.9405;1.9725;
    2.0205;2.0605;2.1085;1.7965;1.5485;1.3245;1.0925;0.8525;0.6205;1.8445;
    1.5725;1.3085;1.0365;0.7645;0.4925;0.5645]

P = [15,0,33;
    12,0,33;
    9,0,33;
    6,0,33;
    3,0,33;
    0,3,33;
    0,6,33;
    0,9,33;
    0,12,33;
    0,15,33;
    0,18,33;
    15,0,36;
    12,0,36;
    9,0,36;
    6,0,36;
    3,0,36;
    0,3,36;
    0,6,36;
    0,9,36;
    0,12,36;
    0,15,36;
    0,18,36;
    15,0,39;
    15,0,42;
    15,0,45;
    15,0,48;
    15,0,51;
    15,0,54;
    0,18,39;
    0,18,42;
    0,18,45;
    0,18,48;
    0,18,51;
    0,18,54]

L=zeros(1,12);
% U=zeros(1,1);

% P=[15,0,36;
%    12,0,36;
%    9,0,36;
%    6,0,36;
%    3,0,36;
%    0,3,36;
%    0,6,36;
%    0,9,36;
%    0,12,36;
%    0,15,36;
%    0,18,36]
% L=zeros(1,12)

% for i = 1:length(x)-1
%     L1=[P(i,1),P(i,2),P(i,3), 0,0,0,-x(i)*P(i,1),-x(i)*P(i,2),-x(i)*P(i,3)]
%     L2= [0,0,0, -y(i)*P(i,1),-y(i)*P(i,2),-y(i)*P(i,3),P(i,1),P(i,2),P(i,3)]
%     LP=[L1;L2]
%     L=[L;LP]
%     U=[U;x(i);y(i)]
% end
for i = 1:length(x)-1
    L1=[P(i,1),P(i,2),P(i,3),1,0,0,0,0, -x(i)*P(i,1),-x(i)*P(i,2),-x(i)*P(i,3),-x(i)]
    L2= [0,0,0,0,P(i,1),P(i,2),P(i,3),1,-y(i)*P(i,1),-y(i)*P(i,2),-y(i)*P(i,3),-y(i)]
    L=[L;L1;L2]
end
[U,S,V]=svd(L)
X = V(:,end)
a1=[X(1);X(2);X(3)]
a2=[X(5);X(6);X(7)]
a3=[X(9);X(10);X(11)]

rho = -1 / norm(a3)
r3 = rho*a3
u_0 = rho^2*(a1'*a3)
v_0 = rho^2*(a2'*a3)

theta = acosd(-(cross(a1,a3)'*cross(a2,a3))/(norm(cross(a1,a3))'*norm(cross(a2,a3))))
alpha = rho^2 * norm(cross(a1,a3))*sind(theta)
beta = rho^2 * norm(cross(a2,a3))*sind(theta)

r1= (cross(a2,a3)/norm(cross(a2,a3)))
r3=a3
r2=cross(r3,r1)

% for i = 1:length(x)-1
%     L1=[P(i,1),P(i,2),P(i,3), -x(i)*P(i,1),-x(i)*P(i,2),-x(i)*P(i,3),0,0,0,1,0]
%     L2= [0,0,0, -y(i)*P(i,1),-y(i)*P(i,2),-y(i)*P(i,3),P(i,1),P(i,2),P(i,3),0,1]
%     LP=[L1;L2]
%     L=[L;LP]
%     U=[U;x(i);y(i)]
% end

% for i = 1:length(x)-1
%     L1=[P(i,1),P(i,2),P(i,3),1,0,0,0,0, -x(i)*P(i,1),-x(i)*P(i,2),-x(i)*P(i,3),-x(i)]
%     L2= [0,0,0,0,P(i,1),P(i,2),P(i,3),1,-y(i)*P(i,1),-y(i)*P(i,2),-y(i)*P(i,3),-y(i)]
%     LP=[L1;L2]
%     L=[L;LP]
%     U=[U;x(i);y(i)]
% end

% [u,s,v]=svd(L)
% X= v
% %X = inv(L'*L)*L'*U
% 
% T1=[X(1);X(2);X(3)]
% T2=[X(4);X(5);X(6)]
% T3=[X(7);X(8);X(9)]
% C1=X(10)
% C2=X(11)
% 
% nt1 = norm(a1)
% nt2 = norm(a2)
% nt3 = norm(a3)
% 
% alpha = norm(cross(T1,T2))/nt2^2
% beta =  norm(cross(T2,T3))/nt2^2
% u0 = (T1'*T2)/nt2^2
% v0 = (T2'*T2)/nt2^2
% 
% tx= (nt2/(norm(cross(T1,T2))))*(C1-((T1'*T2))/nt2^2)
% ty= (nt2/(norm(cross(T2,T3))))*(C2-((T2'*T3))/nt2^2)
% tz=1/nt2
% 
% r1= (nt2/(norm(cross(a1,a2))))*(a1-((a1'*a2)/nt2^2)*a2)
% r2= (nt2/(norm(cross(a2,a3))))*(a3-((a2'*a3)/nt2^2)*a2)
% r3=a2/nt2
% 
rotbeta = asind(-r3(1))
rotalpha = asind(r3(2)/cosd(rotbeta))
rotgamma = atand(r2(1)/r1(1))
k = [alpha, -alpha*(cosd(theta)/sind(theta)),u_0;
     0,beta/sin(theta), v_0;
     0,0,1];
t = rho*inv(k)*[X(4);X(8);X(12)] 

phi_rot = -atand(r2(3)/r3(3))
gamma_rot = -atand(r1(2)/r1(1))
omega_rot = atand(r1(3)/(-r2(3)*sin(phi_rot)+r3(3)*cos(phi_rot)))