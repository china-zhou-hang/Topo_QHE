clear all;
clc;
addpath(genpath('D:\Matlab2021b\DrosteEffect-BrewerMap-3.2.8.0'));
assert(~isempty(which('brewermap')), 'brewermap 未在路径上');

cmap = brewermap(256,'-RdYlBu');   % 需要反向就用 '-RdYlBu'
set(groot,'DefaultFigureColormap',cmap);  % 之后新建的所有 figure 默认用此色表
nh=10000;
rh = 0.0016;
rc = 18;
Gc= 200;
kB=0.6950;
Ts=6000;
Ta=300;
E1 = 0; E2 = 1451; E3 = 13229; E4 = 14840;  
E41=(E4-E1);
Glist   = logspace(log10(1e1), log10(1e6), 100);
% Omelist = logspace(-2, log10(30), 100)*1611;
Omelist = linspace(0, 0.6, 100)*1611;
Vlist = (208.5 * log(200./Glist) + 11778)/11778;
[X, Y] = meshgrid(Omelist/1611, Vlist);

Omelist1 = linspace(0, 0.6, 10)*1611;

a3 = 0.5*Glist;
a4 = 0.5*(rh*(nh+1) + rc);
OmeTline=sqrt(a3.*a4)/1611;

V=0.95*11778;
G = 200*exp(11778*(1-V/11778)/208.5);
OmeT = sqrt(0.5*G*(0.5*(rh*(nh+1) + rc)))/1611;
etaC=0.95;

% load('fig2data/NR_eta0.mat')%1-1
load('fig2data/RWA_eta.mat')%1-1
load('fig2data/RWA_rho43.mat')%1-2

load('fig2data/RWA_eta_1D.mat')%2-1
load('fig2data/NR_eta_1D01.mat')%2-1
load('fig2data/RWA_detadO2_1D.mat')%2-1

load('fig2data/RWA_rho43_1D.mat')%2-3
load('fig2data/RWA_Jcoh_1D.mat')%2-3

etamin = G*Gc*V*rc/(E41*etaC*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh));
etamax = G*Gc*V/(E41*etaC*(G*(Gc+nh*rh)+3*nh*rh*Gc));

% disp(etamin)
% disp(etamax)
% disp(max(2*RWA_rho43(:)))
% disp(min(2*RWA_rho43(:)))
disp(max(RWA_eta(:)))
disp(min(RWA_eta(:)))
disp((max(RWA_eta(:))-min(RWA_eta(:)))/min(RWA_eta(:)))
%%
% clc
close all
brcolor=[237,249,249]/255;
excolor=[254,245,232]/255;
FaceA=1;
textbrcolor=[0 0 0];
textexcolor=[0 0 0];
lft=0.1;
hgap = 0.1;  % 水平间隙
width = 0.213;  % 子图宽度
height = 0.33; % 子图高度
bot_bot=0.15;
top_bot=0.54;

h = figure('Position', [-1700,50,1600,900]);
ax1 = axes('Position', [lft, top_bot, width, height]); 
contourf(X, Y, RWA_eta, 100, 'LineColor', 'none');hold on;
plot(OmeTline,Vlist,'--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
cb1 = colorbar(ax1,'Location','northoutside');
axpos = ax1.Position;               % 主轴位置
cb1.Position(1) = axpos(1);          % 左对齐
cb1.Position(3) = axpos(3)/2;          % 宽度一致
cb1.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax1, 'FontSize', 18);
clim(ax1, [0.13 0.72]); 
cb1.Ticks = [0.13 0.42 0.72];
cb1.TickLabels = {'', '', ''};
ax1.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
text(0.05,0.94,'i', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.52, 1.06, 'Efficiency $\eta$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(-0.02, 1.15, '$12\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.22, 1.15, '$42\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.45, 1.15, '$72\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
plot(ax1, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax1, linspace(0,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
text(-0.2,1.2, '(a)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);


ax2 = axes('Position', [lft+width+hgap, top_bot, width, height]); 
fill([100*ones(size(Vlist)) fliplr(OmeTline)],[Vlist flip(Vlist)],  excolor, 'FaceAlpha', 1, 'EdgeColor', 'none');hold on;
fill([0*ones(size(Vlist)) OmeTline], [flip(Vlist) Vlist], brcolor, 'FaceAlpha', 1, 'EdgeColor', 'none'); hold on;
plot(OmeTline,Vlist,'--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
set(ax2, 'FontSize', 18);
clim(ax2, [-20 110]); 
ax2.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
plot(ax2, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax2, linspace(0.01,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
text(0.05,0.94,'i', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.1, 0.1, '$\nu=0$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28);
text(0.05, 0.25, 'Trivial', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textbrcolor);
text(0.5, 0.5, '$\nu=1$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28);
text(0.4, 0.7, 'Nontrivial', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);
% text(0.65, 0.55,'TPT','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', -50)
text(-0.2,1.2, '(b)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);


ax3 = axes('Position', [lft+2*width+2*hgap, top_bot, width, height]); 
contourf(X, Y, 2*RWA_rho43, 100, 'LineColor', 'none');hold on;
plot(OmeTline,Vlist,'--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); hold on;
cb3 = colorbar(ax3,'Location','northoutside');
axpos = ax3.Position;               % 主轴位置
cb3.Position(1) = axpos(1);          % 左对齐
cb3.Position(3) = axpos(3)/2;          % 宽度一致
cb3.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax3, 'FontSize', 18);
clim(ax3, [0 0.12]); 
cb3.Ticks = [0 0.06 0.12];
cb3.TickLabels = {'', '', ''};
ax3.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
text(0.05,0.94,'i', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.52, 1.06, 'Coherence $\mathcal{C}_{\ell_1}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.0, 1.15, '$0$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.2, 1.15, '$0.06$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.43, 1.15, '$0.12$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
plot(ax3, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax3, linspace(0,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
text(-0.2,1.2, '(c)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);

ax4 = axes('Position', [lft, bot_bot, width, height]); 
hold(ax4,'on');
yyaxis left
patch(ax4, [0 OmeT OmeT 0], [0.0 0.0 1 1]*100,brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
patch(ax4, [OmeT 1 1 OmeT], [0.0 0.0 1 1]*100,excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
plot(ax4, OmeT*ones(1,20), linspace(0,100,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
plot(ax4, linspace(0,0.6,20),99*etamin*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax4, linspace(0,0.6,20),99*etamax*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);
h1=plot(Omelist/1611, 100*RWA_eta_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[100,136,192]/255);
h2=plot(Omelist1/1611, 100*NR_eta_1D01,  'o','LineWidth',1.8,'MarkerSize', 15, 'Color',[16,166,122]/255,'MarkerFaceColor', [100,136,192]/255);
ax4 = gca;
ax4.YColor = [0 0 0];
ylim(ax4, [20 80]);
yticks(ax4, [20 100*etamin 50 98*etamax 80]);
yticklabels(ax4, {'20%','','50%','','80%'});
ylabel('$\eta$','Interpreter','latex', 'FontSize', 28); 
text(-0.25, 0.2,'Efficiency $\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

yyaxis right
h3=plot(Omelist/1611, 1e5*RWA_detadO2_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[238,90,57]/255);
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax4, [-1 3]);
yticks(ax4, [-1 0 3]);
yticklabels(ax4, {'-1','0','3'});
ax = gca;
ax.YColor = [0,0,0]/255;
text(1.085, 0.55,'($10^{-5}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

ax4.Box      = 'on';
ax4.Layer    = 'top';
ax4.LineWidth = 2;
ax4.TickDir  = 'in';
set(ax4, 'FontSize', 18);
xlim(ax4, [0 0.6]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax4, 0.05, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);

text(-0.2, 0.14, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, 0.88, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
lgd=legend([h1, h2,h3], {'RWA','Exact','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd.Position = [0.20, 0.27, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 22;


ax5 = axes('Position', [lft+width+hgap, bot_bot, width, height]); 
hold(ax5,'on');
nu = ones(size(Omelist));
nu(Omelist/1611 < OmeT) = 0;
nu(Omelist/1611 >=OmeT) = 1; 
patch(ax5, [0 OmeT OmeT 0], [-1 -1 2 2],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
patch(ax5, [OmeT 4 4 OmeT], [-1 -1 2 2],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
plot(ax5, OmeT*ones(1,20), linspace(-1,2,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
h2=plot(Omelist/1611,nu, 'Color',[100,136,192]/255,'LineWidth', 5);
ax5.Box      = 'on';
ax5.Layer    = 'top';
ax5.LineWidth = 2;
ax5.TickDir  = 'in';
set(ax5, 'FontSize', 18);
xlim(ax5, [0 0.6]);
ylim(ax5, [-1 2]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax5, 0.05, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
% ylabel('$\nu$','Interpreter','latex', 'FontSize', 26); 
text(-0.2, 0,'Winding number $\nu$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
% text(0.05, 1.08, 'Trivial', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textbrcolor);
% text(0.5, 1.08, 'Nontrivial', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);
set([h1 h2],'HandleVisibility','off');




ax6 = axes('Position', [lft+2*width+2*hgap, bot_bot, width, height]); 
hold(ax6,'on');
yyaxis left
patch(ax6, [0 OmeT OmeT 0], [-1 -1 0.12 0.12],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
patch(ax6, [OmeT 4 4 OmeT], [-1 -1 0.12 0.12],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
plot(ax6, OmeT*ones(1,20), linspace(-0.1,0.12,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
h1=plot(Omelist/1611, 2*RWA_rho43_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[100,136,192]/255);
ylim(ax6, [-0.002 0.04]);
yticks(ax6, [0 0.02 0.04]);
yticklabels(ax6, {'0', '0.02','0.04'});
ax = gca;
ax.YColor = [0 0 0];
text(-0.2, 0.15,'Coherence $\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

yyaxis right
h2=plot(Omelist/1611, RWA_Jcoh_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[238,90,57]/255);
ylim(ax6, [-0.8 16]);
yticks(ax6, [0 8 16]);
yticklabels(ax6, {'0', '','16'});
ax.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

ax6.Box      = 'on';
ax6.Layer    = 'top';
ax6.LineWidth = 2;
ax6.TickDir  = 'in';
set(ax6, 'FontSize', 18);
xlim(ax6, [0 0.6]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax6, 0.05, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
lgd=legend([h1, h2], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd.Position = [0.8, 0.2, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 22;