clear all;
clc;
addpath(genpath('D:\Matlab2021b\DrosteEffect-BrewerMap-3.2.8.0'));
assert(~isempty(which('brewermap')), 'brewermap 未在路径上');

cmap = brewermap(256,'-RdYlBu');   % 需要反向就用 '-RdYlBu'
set(groot,'DefaultFigureColormap',cmap);  % 之后新建的所有 figure 默认用此色表
tic
Vlist = linspace(0.8,1.0,1000);
[P43,PEP] = fun01(Vlist);

%%
close all
color1=[218,243,241]/255;
color2=[150,150,150]/255;
figure('Position', [-1000,200,700,400]);
ax1= subplot(1,1,1);
patch(ax1, [0.9 1.1 1.1 0.9], [-10 -10 100 100],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax1, [1.05 1.0 1.0 1.05], [-10 -10 100 100],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
plot(Vlist,PEP,'LineWidth',6, 'Color',[238,90,57]/255);hold on;
plot(0.95*ones(1,20), linspace(-0.5,0.3857,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(linspace(0.6,0.95,20), 0.3857*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(0.95, 1.27, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [238,90,57]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.0 1.05]);
xticklabels({'0.8','0.9','0.95','1.0','1.05'});
ylim([-2 100])
% yticks([0 3 6]);
% yticklabels({'0','3.0','6.0'});
xlabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
ylabel('$\Omega_{\rm EP}$/$\Omega_{\rm c}$','Interpreter','latex', 'FontSize', 24); 
text(0.45,0.15,'(0.95,1.27)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])


function [P43,PEP] = fun01(Vlist)
rh = 0.005;
rc = 140;
Gc= 200;
nh=60000;
nc=0;
Nc=0;
Eac=1611;
Eab=14856;
n=length(Vlist);
eta=zeros(1,n);
jj=zeros(1,n);
P43=zeros(1,n);
PEP=zeros(1,n);
Geff=(rh*(nh+1) + rc*(nc+1));
for i =1:n
    V=Vlist(i);
    G = 200*exp(11778*(1-V)/208.5);
    OmeT = sqrt(0.5 * (G + rc*nc)*0.5 * (rh*(nh+1) + rc*(nc+1)));
    c3=(G+rc+nh*rh)*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh);
    c4=G*(rc+rh*nh)*(G*Gc+G*nh*rh+3*Gc*nh*rh);
    Omec=sqrt(c3/c4)*OmeT;
    OmeEP=abs(Geff-G)/4;
    % PEP(1,i)=Omec/OmeEP;
    PEP(1,i)=OmeEP/Omec;
end
end
