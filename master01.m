clear;
close all;clc
syms rho44 rho43 rho42 rho41 rho34 rho33 rho32 rho31  rho24 rho23 rho22 rho21 rho14 rho13 rho12 rho11
syms rh rc Gc G phi real
syms E1 E2 E3 E4
syms nc nh Nc real
syms Omeab Omeac Ome12 del Del real
phi=0;
p1=[1,0,0,0]';
p2=[0,1,0,0]';
p3=[0,0,1,0]';
p4=[0,0,0,1]';
S23=p2*p3';
S23d=S23';
S14=p1*p4';
S14d=S14';
S34=p3*p4';
S34d=S34';
S12=p1*p2';
S12d=S12';
S13=p1*p3';
S13d=S13';
rho=[rho11,rho12,rho13,rho14;
    rho21,rho22,rho23,rho24;
    rho31,rho32,rho33,rho34;
    rho41,rho42,rho43,rho44];
E4=0;
E3=1611;
E2=1611+11634;
E1=14856;
rh=0.005;
rc=35;
Gc=300;
nh=60000;
nc=0;
Nc=0;
H=Omeac*p4*(p3')+(Omeac*p4*(p3'))';
rhot=0;
rhot=rhot-1i*(H*rho-rho*H);
rhot=rhot+rh/2*((nh+1)*(2*S14*rho*S14d-S14d*S14*rho-rho*S14d*S14));
rhot=rhot+rh/2*(nh*(2*S14d*rho*S14-S14*S14d*rho-rho*S14*S14d));

rhot=rhot+rc/2*((nc+1)*(2*S34*rho*S34d-S34d*S34*rho-rho*S34d*S34));
rhot=rhot+rc/2*(nc*(2*S34d*rho*S34-S34*S34d*rho-rho*S34*S34d));

rhot=rhot+Gc/2*((Nc+1)*(2*S12*rho*S12d-S12d*S12*rho-rho*S12d*S12));
rhot=rhot+Gc/2*(Nc*(2*S12d*rho*S12-S12*S12d*rho-rho*S12*S12d));

rhot=rhot+G/2*(2*S23*rho*S23d-S23d*S23*rho-rho*S23d*S23);


disp([char(rhot(1,1)) '==0'])
disp([char(rhot(1,2)) '==0'])
disp([char(rhot(1,3)) '==0'])
disp([char(rhot(1,4)) '==0'])


disp([char(rhot(2,1)) '==0'])
disp([char(rhot(2,2)) '==0'])
disp([char(rhot(2,3)) '==0'])
disp([char(rhot(2,4)) '==0'])

disp([char(rhot(3,1)) '==0'])
disp([char(rhot(3,2)) '==0'])
disp([char(rhot(3,3)) '==0'])
disp([char(rhot(3,4)) '==0'])

disp([char(rhot(4,1)) '==0'])
disp([char(rhot(4,2)) '==0'])
disp([char(rhot(4,3)) '==0'])
disp([char(rhot(4,4)) '==0'])
