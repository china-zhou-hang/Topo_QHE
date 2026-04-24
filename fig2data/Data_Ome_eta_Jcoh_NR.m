% parallel_1D_scan.m
clear all;
% close all;
clc;
% 保留图窗
tic

%% ---------- basis & operators ----------
e1=[1;0;0;0]; e2=[0;1;0;0]; e3=[0;0;1;0]; e4=[0;0;0;1];
P1=e1*e1'; P2=e2*e2'; P3=e3*e3'; P4=e4*e4';
S14=e1*e4'; S14d=S14';
S34=e3*e4'; S34d=S34';
S12=e1*e2'; S12d=S12';
S23=e2*e3'; S23d=S23';

%% ---------- static parameters ----------
E1 = 0; E2 = 1451; E3 = 13229; E4 = 14840;
H0 = diag([E1,E2,E3,E4]);

nh=10000; rh = 0.0018; rc=18; Gc=200;
nc=0; Nc=0; 
phi=0; 
omega43=(E4-E3);

p.P1=P1; p.P2=P2; p.P3=P3; p.P4=P4;
p.S14=S14; p.S14d=S14d; p.S34=S34; p.S34d=S34d;
p.S12=S12; p.S12d=S12d; p.S23=S23; p.S23d=S23d;
p.rh=rh; p.rc=rc; p.Gc=Gc; p.nh=nh; p.nc=nc; p.Nc=Nc;
p.phi=phi; p.H0=H0; p.omega43=omega43;
V=0.95;
p.Delta = 0;
p.G     = 200*exp(11778*(1-V)/208.5);

% 初态
rho0 = P1; y0 = rho0(:);

%% ---------- time grid & solver ----------
Omelist  = linspace(0,0.6,10)*1611; 
wdrive = p.omega43 + p.Delta;
T = 2*pi/abs(wdrive);

nPeriods_relax = 300;
nPeriods_avg   = 20;
PtsPerPeriod   = 400;
Tend  = (nPeriods_relax+nPeriods_avg)*T;
Nt    = (nPeriods_relax+nPeriods_avg)*PtsPerPeriod;
tspan = linspace(0, Tend, Nt);

opts  = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',T/PtsPerPeriod);

%% ---------- results (prealloc) ----------
nO = numel(Omelist);
NR_rho22_1D  = zeros(1,nO);
NR_rho33_1D  = zeros(1,nO);
NR_rho43_1D  = zeros(1,nO);
NR_Jcoh_1D   = zeros(1,nO);
NR_eta_1D    = zeros(1,nO);

RWA_rho22_1D  = zeros(1,nO);
RWA_rho33_1D  = zeros(1,nO);
RWA_rho43_1D  = zeros(1,nO);
RWA_Jcoh_1D   = zeros(1,nO);
RWA_eta_1D    = zeros(1,nO);

%% ---------- parallel loop over Omega ----------
for iO = 1:nO
    p_loc = p;
    p_loc.Omega = Omelist(iO);

    % ---------- Non-RWA evolution ----------
    [tNR, yNR] = ode15s(@(t,y) rhs_norwa(t,y,p_loc), tspan, y0, opts);

    selN = tNR >= (Tend - 100*T);
    t_last = tNR(selN);
    y_last = yNR(selN,:);
    NN = size(y_last,1);

    rho22 = zeros(NN,1); rho33 = zeros(NN,1); rho44 = zeros(NN,1);
    rho43 = zeros(NN,1); Jcoh  = zeros(NN,1);

    for k = 1:NN
        R = reshape(y_last(k,:),4,4);
        theta = wdrive*t_last(k) + p_loc.phi;
        rho22(k) = real(R(2,2));
        rho33(k) = real(R(3,3));
        rho44(k) = real(R(4,4));
        rho43(k) = abs(R(4,3));                 % = abs(exp(1i*theta)*R(4,3))

        Jcoh(k)  = 4*p_loc.Omega*cos(theta)*imag(R(4,3));
    end

    Vloc = (208.5*log(200/p_loc.G) + 11778)/11778;

    NR_eta_1D(iO)   = p_loc.G*mean(rho33)*Vloc/(0.95*p_loc.nh*p_loc.rh*(E4-E1)/11778);
    NR_rho22_1D(iO) = mean(rho22);
    NR_rho33_1D(iO) = mean(rho33);
    NR_rho43_1D(iO) = mean(rho43);
    NR_Jcoh_1D(iO)  = mean(Jcoh);

    % ---------- RWA evolution ----------
    [tR, yR] = ode15s(@(t,y) rhs_rwa(t,y,p_loc), tspan, y0, opts);

    selR = tR >= (Tend - 10*T);
    t_last = tR(selR);
    y_last = yR(selR,:);
    NN = size(y_last,1);

    rho22 = zeros(NN,1); rho33 = zeros(NN,1); rho44 = zeros(NN,1);
    rho43 = zeros(NN,1); Jcoh  = zeros(NN,1);

    for k = 1:NN
        R = reshape(y_last(k,:),4,4);
        theta = wdrive*t_last(k) + p_loc.phi;
        rho22(k) = real(R(2,2));
        rho33(k) = real(R(3,3));
        rho44(k) = real(R(4,4));
        rho43(k) = abs(R(4,3));
        Jcoh(k)  = 4*p_loc.Omega*cos(theta)*imag(R(4,3));
    end

    Vloc = (208.5*log(200/p_loc.G) + 11778)/11778;

    RWA_eta_1D(iO)   = p_loc.G*mean(rho33)*Vloc/(0.95*p_loc.nh*p_loc.rh*(E4-E1)/11778);
    RWA_rho22_1D(iO) = mean(rho22);
    RWA_rho33_1D(iO) = mean(rho33);
    RWA_rho43_1D(iO) = mean(rho43);
    RWA_Jcoh_1D(iO)  = mean(Jcoh);

    fprintf('Progress: %d/%d\n', iO, nO);
end

toc
%% 
% close all
x = Omelist/1611;
figure('Color','w','Position',[100 120 1200 1400]);

subplot(2,2,1);
plot(x, NR_eta_1D, '-o', 'DisplayName','Non-RWA'); hold on;
plot(x, RWA_eta_1D,'-s', 'DisplayName','RWA');
legend('Location','best'); xlim([min(x) max(x)]);

% subplot(2,2,2);
% plot(x, NR_Jcoh_1D, '-o', 'DisplayName','Non-RWA'); hold on;
% plot(x, RWA_Jcoh_1D,'-s', 'DisplayName','RWA');
% legend('Location','best'); xlim([min(x) max(x)]);
% 
% subplot(2,2,3);
% plot(x, NR_rho43_1D, '-o', 'DisplayName','Non-RWA'); hold on;
% plot(x, RWA_rho43_1D,'-s', 'DisplayName','RWA');
% legend('Location','best'); xlim([min(x) max(x)]);
% 
% subplot(2,2,4);
% relErr_eta = (NR_eta_1D - RWA_eta_1D)./NR_eta_1D;
% plot(x, relErr_eta, '-o', 'DisplayName','(NR-RWA)/NR');
% legend('Location','best'); xlim([min(x) max(x)]);

NR_eta_1D01=NR_eta_1D;
save('NR_eta_1D01.mat', 'NR_eta_1D01');
% save('NR_Jcoh_1D01.mat', 'NR_Jcoh_1D01');
% save('NR_rho43_1D01.mat', 'NR_rho43_1D01');
% save('NR_reta_1D01.mat', 'NR_reta_1D01');
%% ---------- local functions ----------
function dy = rhs_rwa(t,y,p)
    rho = reshape(y,4,4);
    wdrive = p.omega43 + p.Delta;
    Hdrive = p.Omega*(p.S34+ p.S34d);
    Hlab = Hdrive;
    drho = -1i*(Hlab*rho - rho*Hlab);
    drho = drho + dissipators(rho,p);
    dy = drho(:);
end


function dy = rhs_norwa(t,y,p)
    rho = reshape(y,4,4);
    wdrive = p.omega43 + p.Delta;
    Hdrive = 2*p.Omega*cos(wdrive*t + p.phi) * (p.S34 + p.S34d);
    Hlab   = p.H0 + Hdrive;
    drho = -1i*(Hlab*rho - rho*Hlab);
    drho = drho + dissipators(rho,p);
    dy = drho(:);
end


function L = D(rho,A)
    L = A*rho*A' - 0.5*(A'*A*rho + rho*A'*A);
end

function drho = dissipators(rho,p)
    drho = zeros(4,4);
    drho = drho + (p.rh)*((p.nh+1)*D(rho,p.S14) + p.nh*D(rho,p.S14d));
    drho = drho + (p.rc)*((p.nc+1)*D(rho,p.S34) + p.nc*D(rho,p.S34d));
    drho = drho + (p.Gc)*((p.Nc+1)*D(rho,p.S12) + p.Nc*D(rho,p.S12d));
    drho = drho + (p.G/2)*(2*p.S23*rho*p.S23d - p.S23d*p.S23*rho - rho*p.S23d*p.S23);
end
