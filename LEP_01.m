clear all;
% close all;
clc;
syms rhoaa rhoac rhoav rhoab rhoca rhocc rhocv rhocb  rhova rhovc rhovv rhovb rhoba rhobc rhobv rhobb
syms Omeac G positive real
syms Del real
syms nh rh rc Gc positive real
pa=[1,0,0,0]';
pc=[0,1,0,0]';
pv=[0,0,1,0]';
pb=[0,0,0,1]';
Svc=pv*pc';
Svcd=Svc';
Sba=pb*pa';
Sbad=Sba';
Sca=pc*pa';
Scad=Sca';
Sbv=pb*pv';
Sbvd=Sbv';
Sbc=pb*pc';
Sbcd=Sbc';
rho=[rhoaa,rhoac,rhoav,rhoab;
    rhoca,rhocc,rhocv,rhocb;
    rhova,rhovc,rhovv,rhovb;
    rhoba,rhobc,rhobv,rhobb];
nc=0;
Nc=0;
nh=60000;
rh = 0.005;
rc = 140;
Gc= 200;
kB=0.6950;
Ts=6000;
Ta=300;
% Omeac=1e10;
% G=4e2;
%% 

Eac=1611;
Evb=1611;
Ecv=11634;
Eb=0;
Ec=1611;
Ev=1611+11634;
Ea=1611+1611+11634;
H=0;

% H=H+Omeac*pa*(pc')+(Omeac*pa*(pc'))'+Omeac*pv*(pb')+(Omeac*pv*(pb'))';
H=H+Omeac*pa*(pc')+(Omeac*pa*(pc'))';
rhot=0;
rhot=rhot-1i*(H*rho-rho*H);
rhot=rhot+rh/2*((nh+1)*(2*Sba*rho*Sbad-Sbad*Sba*rho-rho*Sbad*Sba));
rhot=rhot+rh/2*(nh*(2*Sbad*rho*Sba-Sba*Sbad*rho-rho*Sba*Sbad));

rhot=rhot+rc/2*((nc+1)*(2*Sca*rho*Scad-Scad*Sca*rho-rho*Scad*Sca));
rhot=rhot+rc/2*(nc*(2*Scad*rho*Sca-Sca*Scad*rho-rho*Sca*Scad));

rhot=rhot+Gc/2*((Nc+1)*(2*Sbv*rho*Sbvd-Sbvd*Sbv*rho-rho*Sbvd*Sbv));
rhot=rhot+Gc/2*(Nc*(2*Sbvd*rho*Sbv-Sbv*Sbvd*rho-rho*Sbv*Sbvd));

rhot=rhot+G/2*(2*Svc*rho*Svcd-Svcd*Svc*rho-rho*Svcd*Svc);

LL=0;
LL=LL-1i*(kron(eye(4),H)-kron(H.',eye(4)));
LL=LL+rh/2*((nh+1)*(2*kron(Sba,conj(Sba))-kron(Sbad*Sba,eye(4))-kron(eye(4),Sba.'*conj(Sba))));
LL=LL+rh/2*(nh*(2*kron(Sbad,conj(Sbad))-kron(Sba*Sbad,eye(4))-kron(eye(4),Sbad.'*conj(Sbad))));
LL=LL+rc/2*((nc+1)*(2*kron(Sca,conj(Sca))-kron(Scad*Sca,eye(4))-kron(eye(4),Sca.'*conj(Sca))));
LL=LL+rc/2*(nc*(2*kron(Scad,conj(Scad))-kron(Sca*Scad,eye(4))-kron(eye(4),Scad.'*conj(Scad))));
LL=LL+Gc/2*((Nc+1)*(2*kron(Sbv,conj(Sbv))-kron(Sbvd*Sbv,eye(4))-kron(eye(4),Sbv.'*conj(Sbv))));
LL=LL+Gc/2*(Nc*(2*kron(Sbvd,conj(Sbvd))-kron(Sbv*Sbvd,eye(4))-kron(eye(4),Sbvd.'*conj(Sbvd))));
LL=LL+G/2*(1*2*kron(Svc,conj(Svc))-kron(Svcd*Svc,eye(4))-kron(eye(4),Svc.'*conj(Svc)));
% rho1=rho(:);
% disp(simplify(rhot-reshape(LL*rho1,4,4)))
disp(LL)
% disp(LL(1:8,1:8))
% disp(LL(9:16,9:16))
% disp(LL(9:16,1:8))
% disp(LL(1:8,9:16))

% rows = [3,9,7,10,12,13,14];
% rows = [2,5,8,14,3,9,12,15,1];
% LL = LL(rows, rows);
% disp(LL)
% AA=eig(LL);
% [V,D]=eig(LL);
%% 
% clc
% aa=diag(D);
% disp(aa(1:11))
% disp(reshape(V(:,1),4,4))
% disp(trace(reshape(V(:,1),4,4)))
% disp(reshape(V(:,1),4,4)/trace(reshape(V(:,1),4,4)))
% % 
% disp(simplify(AA(1)))
% disp(simplify(AA(2)))
% disp(simplify(AA(3)))
% disp(simplify(AA(4)))
% disp(simplify(AA(5)))
% disp(simplify(AA(6)))
% disp(simplify(AA(7)))
% disp(simplify(AA(8)))
% disp(simplify(AA(9)))
% disp(simplify(AA(10)))
% disp(simplify(AA(11)))
% disp(simplify(AA(12)))
% disp(simplify(AA(13)))
% disp(simplify(AA(14)))
% disp(simplify(AA(15)))
% disp(simplify(AA(16)))
% % 
% % disp(AA(1))
% % disp(AA(2))
% % disp(AA(3))
% % disp(AA(4))
% % disp(AA(5))
% % disp(AA(6))
% % disp(AA(7))
% % disp(AA(8))
% % disp(AA(9))
% % disp(AA(10))
% % disp(AA(11))
% % disp(AA(12))
% % disp(AA(13))
% % disp(AA(14))
% % disp(AA(15))
% % disp(AA(16))
