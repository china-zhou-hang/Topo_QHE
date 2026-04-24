% parallel_1D_scan.m
clear all;close all;clc;
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

nh=10000; rh = 0.01; rc=100; Gc=200;
nc=0; Nc=0; 
phi=0; 
omega43=(E4-E3);

p.P1=P1; p.P2=P2; p.P3=P3; p.P4=P4;
p.S14=S14; p.S14d=S14d; p.S34=S34; p.S34d=S34d;
p.S12=S12; p.S12d=S12d; p.S23=S23; p.S23d=S23d;
p.rh=rh; p.rc=rc; p.Gc=Gc; p.nh=nh; p.nc=nc; p.Nc=Nc;
p.phi=phi; p.H0=H0; p.omega43=omega43;
V=0.95;
G     = 200*exp(11778*(1-V)/208.5);

Omelist  = linspace(0,0.6,100)*1611;
% Omelist  = logspace(-2, log10(30),100)*1611;
%% ---------- results (prealloc) ----------
nO = numel(Omelist);
RWA_rho22_1D = zeros(1,nO);
RWA_rho33_1D = zeros(1,nO);
RWA_rho43_1D = zeros(1,nO);
RWA_Jcoh_1D  = zeros(1,nO);
RWA_eta_1D   = zeros(1,nO);
RWA_Psun_1D  = zeros(1,nO);
RWA_detadO2_1D  = zeros(1,nO);
for iO = 1:nO
    Omeac = Omelist(iO);
    rho44=(1.0e+4*(625.0*G^2 + 2.13e+4*G + 2500.0*Omeac^2))/(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9);
    rho33=(2.0e+4*(1250.0*Omeac^2 + 5620.0*G + 1.91e+5))/(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9);
    rho43=(1.25e+7*(18.0*Omeac - 1.0*G*Omeac))/(G^2*2.01e+7i + G*Omeac^2*1.69e+6i + G*7.96e+8i + Omeac^2*7.5e+7i + 3.83e+9i);
    RWA_eta_1D(iO) = (0.088569566079981860952866819714877*G*((417*log(200/G))/2 + 11778)*(1250.0*Omeac^2 + 5625.0*G + 191259.0))/(20094375.0*G^2 + 1687500.0*G*Omeac^2 + 795740901.0*G + 75002500.0*Omeac^2 + 3825180000.0);
    RWA_rho43_1D(iO)=abs(rho43);
    RWA_Jcoh_1D(iO) =2*Omeac*abs(rho43);
    RWA_detadO2_1D(iO)=(7812500*G*((417*log(200/G))/2 + 11778))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)) - (15625000*G*Omeac*(150005000*Omeac + 3375000*G*Omeac)*((417*log(200/G))/2 + 11778))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)^2) + (6250*G*(150005000*Omeac + 3375000*G*Omeac)^2*((417*log(200/G))/2 + 11778)*(1250*Omeac^2 + 5625*G + 191259))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)^3) - (3125*G*(3375000*G + 150005000)*((417*log(200/G))/2 + 11778)*(1250*Omeac^2 + 5625*G + 191259))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)^2);
end

toc
%% 
close all
x = Omelist/1611;  % 无量纲化横轴
% x = 10.^(-2 + 2*x); 

figure('Color','w','Position',[-1600 120 1200 500]);

ax1=subplot(1,2,1);
plot(x, RWA_eta_1D,  '-s', 'DisplayName','RWA'); grid on; box on;
legend('Location','best'); xlim([min(x) max(x)]);
% set(gca, 'XScale', 'log');
% xlim(ax1, [0.01 30]);
% xticks([0.01 0.1 0.3780 1 5 30]);
% xticklabels({'','0.1','','1','10','30'});
% text(-0.01, -0.06, '0^{+}', 'Units', 'normalized', 'FontSize', 18);
% xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
% ylabel('$\eta$','Interpreter','latex', 'FontSize', 26); 
% text(0.1, 0.8,'V=0.8','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)

ax2=subplot(1,2,2);
yyaxis left
plot(x, RWA_rho43_1D,  '-s', 'DisplayName','RWA'); grid on; box on;
% legend('Location','best'); xlim([min(x) max(x)]);
% set(gca, 'XScale', 'log');
% ylabel('$\rho_{43}$','Interpreter','latex', 'FontSize', 26); 

yyaxis right
plot(x, RWA_Jcoh_1D,  '-s', 'DisplayName','RWA'); grid on; box on;
% legend('Location','best'); xlim([min(x) max(x)]);
% set(gca, 'XScale', 'log');
% xlim(ax2, [0.01 30]);
% xticks([0.01 0.1 0.3780 1 5 30]);
% xticklabels({'','0.1','','1','10','30'});
% text(-0.01, -0.06, '0^{+}', 'Units', 'normalized', 'FontSize', 18);
% ylabel('$J_c$','Interpreter','latex', 'FontSize', 26); 

xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
save('RWA_eta_1D.mat', 'RWA_eta_1D');
save('RWA_Jcoh_1D.mat', 'RWA_Jcoh_1D');
save('RWA_rho43_1D.mat', 'RWA_rho43_1D');
save('RWA_detadO2_1D.mat', 'RWA_detadO2_1D');
% disp(max(RWA_eta_1D(:)))
% disp(min(RWA_eta_1D(:)))