clear;
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
nh=10000;
rh = 0.0016;
rc = 18;
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
%%
% [V,D]=eig(LL);
rho_ss_vec = null(LL);
N = sqrt(size(LL,1));
rho_ss = reshape(rho_ss_vec, N, N);
rho_ss = rho_ss / trace(rho_ss);
%%
clc
% disp(simplify(rho_ss))
Von=-trace(rho_ss*logm(rho_ss));
disp(vpa(simplify(Von),3))
jj=G*rho_ss(2,2);
VV=Ecv+kB*Ta*log(rho_ss(2,2)/rho_ss(3,3));
PP=jj*VV;
disp(jj)
disp(VV)
disp(PP)