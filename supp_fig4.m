clear all;
clc;
addpath(genpath('D:\Matlab2021b\DrosteEffect-BrewerMap-3.2.8.0'));
assert(~isempty(which('brewermap')), 'brewermap 未在路径上');

cmap = brewermap(256,'-RdYlBu');   % 需要反向就用 '-RdYlBu'
set(groot,'DefaultFigureColormap',cmap);  % 之后新建的所有 figure 默认用此色表
nh=10000;
rh = 0.0018;
rc = 18;
Gc= 200;
kB=0.6950;
Ts=6000;
Ta=300;
E1 = 0; E2 = 1451; E3 = 13229; E4 = 14840;  
E41=(E4-E1);
V1=0.95*11778;
G1 = 200*exp(11778*(1-V1/11778)/208.5);

OmeT = sqrt(0.5*G1*(0.5*(rh*(nh+1) + rc)))/11778;
Omestar=sqrt(OmeT^2+((G1+rh*(nh+1) + rc)/4)^2)/11778;
load('supp_fig4_time_plot/tlist.mat')%1-1

load('supp_fig4_time_plot/p44_1.mat')%1-1
load('supp_fig4_time_plot/p44_2.mat')%1-1
load('supp_fig4_time_plot/p44_3.mat')%1-1
load('supp_fig4_time_plot/p44_4.mat')%1-1
load('supp_fig4_time_plot/p44_5.mat')%1-1

load('supp_fig4_time_plot/p33_1.mat')%1-1
load('supp_fig4_time_plot/p33_2.mat')%1-1
load('supp_fig4_time_plot/p33_3.mat')%1-1
load('supp_fig4_time_plot/p33_4.mat')%1-1
load('supp_fig4_time_plot/p33_5.mat')%1-1

load('supp_fig4_time_plot/p22_1.mat')%1-1
load('supp_fig4_time_plot/p22_2.mat')%1-1
load('supp_fig4_time_plot/p22_3.mat')%1-1
load('supp_fig4_time_plot/p22_4.mat')%1-1
load('supp_fig4_time_plot/p22_5.mat')%1-1

load('supp_fig4_time_plot/p11_1.mat')%1-1
load('supp_fig4_time_plot/p11_2.mat')%1-1
load('supp_fig4_time_plot/p11_3.mat')%1-1
load('supp_fig4_time_plot/p11_4.mat')%1-1
load('supp_fig4_time_plot/p11_5.mat')%1-1

load('supp_fig4_time_plot/p4p3_1.mat')%1-1
load('supp_fig4_time_plot/p4p3_2.mat')%1-1
load('supp_fig4_time_plot/p4p3_3.mat')%1-1
load('supp_fig4_time_plot/p4p3_4.mat')%1-1
load('supp_fig4_time_plot/p4p3_5.mat')%1-1

load('supp_fig4_time_plot/p43_1.mat')%1-1
load('supp_fig4_time_plot/p43_2.mat')%1-1
load('supp_fig4_time_plot/p43_3.mat')%1-1
load('supp_fig4_time_plot/p43_4.mat')%1-1
load('supp_fig4_time_plot/p43_5.mat')%1-1
pp1=zeros(length(tlist),6);
pp1(:,1)=tlist';
pp1(:,2)=p44_1';
pp1(:,3)=p44_2';
pp1(:,4)=p44_3';
pp1(:,5)=p44_4';
pp1(:,6)=p44_5';

pp2=zeros(length(tlist),6);
pp2(:,1)=tlist';
pp2(:,2)=p33_1';
pp2(:,3)=p33_2';
pp2(:,4)=p33_3';
pp2(:,5)=p33_4';
pp2(:,6)=p33_5';


pp3=zeros(length(tlist),6);
pp3(:,1)=tlist';
pp3(:,2)=p4p3_1';
pp3(:,3)=p4p3_2';
pp3(:,4)=p4p3_3';
pp3(:,5)=p4p3_4';
pp3(:,6)=p4p3_5';


pp4=zeros(length(tlist),6);
pp4(:,1)=tlist';
pp4(:,2)=p22_1';
pp4(:,3)=p22_2';
pp4(:,4)=p22_3';
pp4(:,5)=p22_4';
pp4(:,6)=p22_5';

pp5=zeros(length(tlist),6);
pp5(:,1)=tlist';
pp5(:,2)=p11_1';
pp5(:,3)=p11_2';
pp5(:,4)=p11_3';
pp5(:,5)=p11_4';
pp5(:,6)=p11_5';


pp6=zeros(length(tlist),6);
pp6(:,1)=tlist';
pp6(:,2)=p43_1';
pp6(:,3)=p43_2';
pp6(:,4)=p43_3';
pp6(:,5)=p43_4';
pp6(:,6)=p43_5';


%% 
close all
color0=[0,0,0]/255;
color1=[255,0,0]/255;
color2=[100,136,192]/255;
color3=[0,176,80]/255;
color4=[119,45,133]/255;
lft=0.1;
hgap = 0.1;  % 水平间隙
width = 0.213;  % 子图宽度
height = 0.32; % 子图高度
bot_bot=0.15;
top_bot=0.6;

h = figure('Position', [-1700,50,1600,900]);
ax1 = axes('Position', [lft, top_bot, width, height]); 
h1=plot(pp1(:,1),pp1(:,2),'LineWidth',6, 'Color', color0);hold on;
h2=plot(pp1(:,1),pp1(:,3),'LineWidth',6, 'Color', color1);hold on;
h3=plot(pp1(:,1),pp1(:,4),'LineWidth',6, 'Color', color2);hold on;
h4=plot(pp1(:,1),pp1(:,5),'LineWidth',6, 'Color', color3);hold on;
h5=plot(pp1(:,1),pp1(:,6),'LineWidth',6, 'Color', color4);hold on;
set(gca, 'XScale', 'log');
xlim(ax1, [0.01 4000]);
ylim(ax1, [-0.02 0.4]);
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 18);
xticks([0.01 1 100 4000]);
xticklabels({'0','1','100','4000'});
yticks(ax1, [0 0.1 0.2 0.3 0.4]);
text(ax1, 0.05, 0.92, '(a)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('Time (fs)','Interpreter','latex', 'FontSize', 28);
ylabel('$\rho_{33}$','Interpreter','latex', 'FontSize', 28); 

ax2 = axes('Position', [lft+width+hgap, top_bot, width, height]); 
h1=plot(pp2(:,1),pp2(:,2),'LineWidth',6, 'Color', color0);hold on;
h2=plot(pp2(:,1),pp2(:,3),'LineWidth',6, 'Color', color1);hold on;
h3=plot(pp2(:,1),pp2(:,4),'LineWidth',6, 'Color', color2);hold on;
h4=plot(pp2(:,1),pp2(:,5),'LineWidth',6, 'Color', color3);hold on;
h5=plot(pp2(:,1),pp2(:,6),'LineWidth',6, 'Color', color4);hold on;
set(gca, 'XScale', 'log');
xlim(ax2, [0.01 4000]);
ylim(ax2, [-0.0003 0.006]);
ax2.Box      = 'on';
ax2.Layer    = 'top';
ax2.LineWidth = 2;
ax2.TickDir  = 'in';
set(ax2, 'FontSize', 18);
xticks([0.01 1 10 100 1000]);
xticklabels({'0','1','100','4000'});
yticks(ax2, [0 0.003 0.006]);
text(ax2, 0.05, 0.92, '(b)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('Time (fs)','Interpreter','latex', 'FontSize', 28);
ylabel('$\rho_{22}$','Interpreter','latex', 'FontSize', 28); 


ax3 = axes('Position', [lft+2*width+2*hgap, top_bot, width, height]); 
h1=plot(pp3(:,1),pp3(:,2),'LineWidth',6, 'Color', color0);hold on;
h2=plot(pp3(:,1),pp3(:,3),'LineWidth',6, 'Color', color1);hold on;
h3=plot(pp3(:,1),pp3(:,4),'LineWidth',6, 'Color', color2);hold on;
h4=plot(pp3(:,1),pp3(:,5),'LineWidth',6, 'Color', color3);hold on;
h5=plot(pp3(:,1),pp3(:,6),'LineWidth',6, 'Color', color4);hold on;
set(gca, 'XScale', 'log');
xlim(ax3, [0.01 4000]);
ylim(ax3, [-0.02 0.4]);
ax3.Box      = 'on';
ax3.Layer    = 'top';
ax3.LineWidth = 2;
ax3.TickDir  = 'in';
set(ax3, 'FontSize', 18);
xticks([0.01 1 100 4000]);
xticklabels({'0','1','100','4000'});
yticks(ax3, [0 0.1 0.2 0.3 0.4]);
yticklabels({'0','0.1','0.2','0.3','0.4'});
text(ax3, 0.05, 0.92, '(c)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('Time (fs)','Interpreter','latex', 'FontSize', 28);
ylabel('$\rho_{33}$-$\rho_{22}$','Interpreter','latex', 'FontSize', 28); 

ax4 = axes('Position', [lft, bot_bot, width, height]); 
h1=plot(pp4(:,1),pp4(:,2),'LineWidth',6, 'Color', color0);hold on;
h2=plot(pp4(:,1),pp4(:,3),'LineWidth',6, 'Color', color1);hold on;
h3=plot(pp4(:,1),pp4(:,4),'LineWidth',6, 'Color', color2);hold on;
h4=plot(pp4(:,1),pp4(:,5),'LineWidth',6, 'Color', color3);hold on;
h5=plot(pp4(:,1),pp4(:,6),'LineWidth',6, 'Color', color4);hold on;
set(gca, 'XScale', 'log');
xlim(ax4, [0.01 4000]);
ylim(ax4, [-0.005 0.1]);
ax4.Box      = 'on';
ax4.Layer    = 'top';
ax4.LineWidth = 2;
ax4.TickDir  = 'in';
set(ax4, 'FontSize', 18);
xticks([0.01 1 100 4000]);
xticklabels({'0','1','100','4000'});
yticks(ax4, [0 0.05 0.1]);
yticklabels({'0','0.05','0.1'});
text(ax4, 0.05, 0.92, '(d)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('Time (fs)','Interpreter','latex', 'FontSize', 28);
ylabel('$\rho_{11}$','Interpreter','latex', 'FontSize', 28); 

ax5 = axes('Position', [lft+width+hgap, bot_bot, width, height]); 
h1=plot(pp5(:,1),pp5(:,2),'LineWidth',6, 'Color', color0);hold on;
h2=plot(pp5(:,1),pp5(:,3),'LineWidth',6, 'Color', color1);hold on;
h3=plot(pp5(:,1),pp5(:,4),'LineWidth',6, 'Color', color2);hold on;
h4=plot(pp5(:,1),pp5(:,5),'LineWidth',6, 'Color', color3);hold on;
h5=plot(pp5(:,1),pp5(:,6),'LineWidth',6, 'Color', color4);hold on;
set(gca, 'XScale', 'log');
xlim(ax5, [0.01 4000]);
ylim(ax5, [0.6 1.02]);
ax5.Box      = 'on';
ax5.Layer    = 'top';
ax5.LineWidth = 2;
ax5.TickDir  = 'in';
set(ax5, 'FontSize', 18);
xticks([0.01 1 100 4000]);
xticklabels({'0','1','100','4000'});
yticks(ax5, [0.6 0.8 1.0]);
yticklabels({'0.6','0.8','1.0'});
text(ax5, 0.05, 0.92, '(e)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('Time (fs)','Interpreter','latex', 'FontSize', 28);
ylabel('$\rho_{00}$','Interpreter','latex', 'FontSize', 28); 

ax6 = axes('Position', [lft+2*width+2*hgap, bot_bot, width, height]); 
h1=plot(pp6(:,1),2*pp6(:,2),'LineWidth',6, 'Color', color0);hold on;
h2=plot(pp6(:,1),2*pp6(:,3),'LineWidth',6, 'Color', color1);hold on;
h3=plot(pp6(:,1),2*pp6(:,4),'LineWidth',6, 'Color', color2);hold on;
h4=plot(pp6(:,1),2*pp6(:,5),'LineWidth',6, 'Color', color3);hold on;
h5=plot(pp6(:,1),2*pp6(:,6),'LineWidth',6, 'Color', color4);hold on;
set(gca, 'XScale', 'log');
xlim(ax6, [0.01 4000]);
ylim(ax6, [-0.005 0.10]);
ax6.Box      = 'on';
ax6.Layer    = 'top';
ax6.LineWidth = 2;
ax6.TickDir  = 'in';
set(ax6, 'FontSize', 18);
xticks([0.01 1 100 4000]);
xticklabels({'0','1','100','4000'});
yticks(ax6, [0 0.05 0.1]);
text(ax6, 0.05, 0.92, '(f)', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('Time (fs)','Interpreter','latex', 'FontSize', 28);
ylabel('$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'FontSize', 28); 

lgd=legend([h1, h2,h3,h4,h5], {'$\Omega=0$','$\Omega=0.5\Omega_{T}$','$\Omega=\Omega_{\rm T}$','$\Omega=5\Omega_{\rm T}$','$\Omega=10\Omega_{\rm T}$'}, 'Interpreter', 'latex');
lgd.Position = [0.11, 0.72, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 18;
