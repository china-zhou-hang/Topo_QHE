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

load('figs6data/RWA_PP.mat')%1-1
load('figs6data/RWA_SvN.mat')%1-2


load('figs6data/dPPdO2_1D.mat')%2-1
load('figs6data/dSvNdO1_1D.mat')%2-1
load('figs6data/PP_1D.mat')%2-1
load('figs6data/SvN_1D.mat')%2-3


disp(max(RWA_PP(:)))
disp(min(RWA_PP(:)))

disp(max(RWA_SvN(:)))
disp(min(RWA_SvN(:)))
%%
close all
brcolor=[237,249,249]/255;
excolor=[254,245,232]/255;
lft=0.11;
hgap = 0.15;  % 水平间隙
width = 0.34;  % 子图宽度
height = 0.31; % 子图高度
bot_bot=0.15;
top_bot=0.54;
h = figure;
set(h, 'Position', [-1200, 0, 1000*0.9, 1000*0.9]); % [left, bottom, width, height]
ax1 = axes('Position', [lft, top_bot, width, height]);
contourf(X, Y, RWA_PP, 100, 'LineColor', 'none');hold on;
plot(OmeTline,Vlist,'--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
cb1 = colorbar(ax1,'Location','northoutside');
axpos = ax1.Position;               % 主轴位置
cb1.Position(1) = axpos(1);          % 左对齐
cb1.Position(3) = axpos(3)/2;          % 宽度一致
cb1.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax1, 'FontSize', 18);
clim(ax1, [4 14]); 
cb1.Ticks = [4 9 14];
cb1.TickLabels = {'', '', ''};
ax1.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
text(0.05,0.94,'(a)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.52, 1.06, 'Power $P/\omega_{32}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(-0.02, 1.15, '$4$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.22, 1.15, '$9$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.45, 1.15, '$14$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
plot(ax1, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax1, linspace(0,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
% text(-0.2,1.2, '(a)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=[0,0,0]);


ax2 = axes('Position', [lft+width+hgap, top_bot, width, height]);
contourf(X, Y, RWA_SvN, 100, 'LineColor', 'none');hold on;
plot(OmeTline,Vlist,'--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
cb2 = colorbar(ax2,'Location','northoutside');
axpos = ax2.Position;               % 主轴位置
cb2.Position(1) = axpos(1);          % 左对齐
cb2.Position(3) = axpos(3)/2;          % 宽度一致
cb2.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax2, 'FontSize', 18);
clim(ax2, [0.3 1.15]); 
cb2.Ticks = [0.3 0.7 1.1];
cb2.TickLabels = {'', '', ''};
ax2.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
text(0.05,0.94,'(b)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.52, 1.06, 'Entropy $S_{\rm vN}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(-0.02, 1.15, '$0.3$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.22, 1.15, '$0.7$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.45, 1.15, '$1.1$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
plot(ax2, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax2, linspace(0,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
% text(-0.2,1.2, '(b)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=[0,0,0]);

ax3 = axes('Position', [lft, bot_bot, width, height]);
hold(ax3,'on');
yyaxis left
patch(ax3, [0 OmeT OmeT 0], [0.0 0.0 1 1]*100,brcolor, 'FaceAlpha',1, 'EdgeColor','none');
patch(ax3, [OmeT 1 1 OmeT], [0.0 0.0 1 1]*100,excolor, 'FaceAlpha',1, 'EdgeColor','none');
plot(ax3, OmeT*ones(1,20), linspace(0,100,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
h1=plot(Omelist/1611, PP_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[100,136,192]/255);
ax3 = gca;
ax3.YColor = [0 0 0];
ylim(ax3, [4 14]);
yticks(ax3, [4 9 14]);
yticklabels(ax3, {'4','9','14'});
% ylabel('$\eta$','Interpreter','latex', 'FontSize', 28); 
text(-0.2, 0.2,'Power $P/\omega_{32}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

yyaxis right
h3=plot(Omelist/1611, dPPdO2_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[238,90,57]/255);
text(1.08, 0.4,'$\mathcal{H}_{P}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax3, [-2 6]);
yticks(ax3, [-2 0 6]);
yticklabels(ax3, {'-2','0','6'});
ax = gca;
ax.YColor = [0,0,0]/255;

ax3.Box      = 'on';
ax3.Layer    = 'top';
ax3.LineWidth = 2;
ax3.TickDir  = 'in';
set(ax3, 'FontSize', 18);
xlim(ax3, [0 0.6]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax3, 0.05, 0.94, '(c)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);

text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
% lgd=legend([h1, h2,h3], {'RWA','Exact','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
% lgd.Position = [0.20, 0.27, 0.1, 0.1];  
% lgd.Box = 'off';
% lgd.FontSize = 22;

ax4 = axes('Position', [lft+width+hgap, bot_bot, width, height]);
hold(ax4,'on');
yyaxis left
patch(ax4, [0 OmeT OmeT 0], [0.0 0.0 1 1]*100,brcolor, 'FaceAlpha',1, 'EdgeColor','none');
patch(ax4, [OmeT 1 1 OmeT], [0.0 0.0 1 1]*100,excolor, 'FaceAlpha',1, 'EdgeColor','none');
plot(ax4, OmeT*ones(1,20), linspace(0,100,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
h1=plot(Omelist/1611, SvN_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[100,136,192]/255);
ax4 = gca;
ax4.YColor = [0 0 0];
ylim(ax4, [0.3 0.9]);
yticks(ax4, [0.3 0.6 0.9]);
yticklabels(ax4, {'0.3','0.6','0.9'});
% ylabel('$\eta$','Interpreter','latex', 'FontSize', 28); 
text(-0.2, 0.2,'Entropy $S_{\rm vN}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

yyaxis right
h3=plot(Omelist/1611, 1e3*dSvNdO1_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[238,90,57]/255);
text(1.08, 0.12,'$\partial S_{vN}/\partial \Omega$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax4, [-1 0]);
yticks(ax4, [-1 0]);
yticklabels(ax4, {'-1','0'});
ax = gca;
ax.YColor = [0,0,0]/255;
text(1.085, 0.65,'($10^{-3}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

ax4.Box      = 'on';
ax4.Layer    = 'top';
ax4.LineWidth = 2;
ax4.TickDir  = 'in';
set(ax4, 'FontSize', 18);
xlim(ax4, [0 0.6]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax4, 0.05, 0.94, '(d)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);

text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
% lgd=legend([h1, h2,h3], {'RWA','Exact','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
% lgd.Position = [0.20, 0.27, 0.1, 0.1];  
% lgd.Box = 'off';
% lgd.FontSize = 22;