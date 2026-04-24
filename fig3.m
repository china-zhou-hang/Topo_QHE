clear all;
clc;
addpath(genpath('D:\Matlab2021b\DrosteEffect-BrewerMap-3.2.8.0'));
assert(~isempty(which('brewermap')), 'brewermap 未在路径上');

cmap = brewermap(256,'-RdYlBu');   % 需要反向就用 '-RdYlBu'
set(groot,'DefaultFigureColormap',cmap);  % 之后新建的所有 figure 默认用此色表
tic
rh = 0.0018;
rc = 18;
Gc= 200;
nh=10000;
Vlist = linspace(0.8,1.0,100);
[P43,PEP,Pstar] = fun01(Vlist,nh,rh,rc,Gc);
%%
close all
color1=[218,243,241]/255;
color2=[150,150,150]/255;
h = figure('Position', [-1000,200,700,600]);
ax1 = axes('Position', [0.2,0.58,0.75,0.35]);
patch(ax1, [0.9 1.1 1.1 0.9], [-0.5 -0.5 6 6],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax1, [1.05 1.0 1.0 1.05], [-0.5 -0.5 6 6],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
plot(Vlist,P43,'LineWidth',6, 'Color',[100,136,192]/255);hold on;
plot(0.95*ones(1,20), linspace(-0.5,0.07,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(linspace(0.6,0.95,20), 0.3857*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.0 1.05]);
xticklabels({''});
ylim([-0.1615 -0.14])
yticks([-0.16 -0.15 -0.14]);
yticklabels({'-0.16','-0.15','-0.14'});
ylabel('$R_{T}$','Interpreter','latex', 'FontSize', 24); 
text(0.9,0.9,'(a)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
plot(0.95, -0.159, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);


ax2 = axes('Position', [0.2,0.15,0.75,0.35]);
patch(ax2, [0.9 1.1 1.1 0.9], [-0.5 -0.5 6 6],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax2, [1.05 1.0 1.0 1.05], [-0.5 -0.5 6 6],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
plot(Vlist,Pstar,'LineWidth',6, 'Color',[238,90,57]/255);hold on;
plot(0.95*ones(1,20), linspace(-0.5,0.07,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(linspace(0.6,0.95,20), 0.3857*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(0.95,0.002, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [238,90,57]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
ax2.Box      = 'on';
ax2.Layer    = 'top';
ax2.LineWidth = 2;
ax2.TickDir  = 'in';
set(ax2, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.0 1.05]);
xticklabels({'0.8','0.9','0.95','1.0','1.05'});
xlabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
ylim([-0.003 0.06])
yticks([0 0.03 0.06]);
yticklabels({'0','0.03','0.06'});
ylabel('$R^{\ast}_{T}$','Interpreter','latex', 'FontSize', 24); 
text(0.9,0.9,'(b)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
% grid on;


function [P43,PEP,Pstar] = fun01(Vlist,nh,rh,rc,Gc)
nc=0;
n=length(Vlist);
P43=zeros(1,n);
PEP=zeros(1,n);
Pstar=zeros(1,n);
Geff=(rh*(nh+1) + rc*(nc+1));
Wh=nh*rh;
for i =1:n
    V=Vlist(i);
    G = 200*exp(11778*(1-V)/208.5);
    Geff=(rh*(nh+1) + rc);
    OmeT = sqrt(0.5 * (G + rc*nc)*0.5 * (rh*(nh+1) + rc*(nc+1)));
    PEP(1,i)=P43(1,i)*OmeT/(abs(Geff-G)/4);

    c3=(G+rc+nh*rh)*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh);
    c4=G*(rc+rh*nh)*(G*Gc+G*nh*rh+3*Gc*nh*rh);
    OmeT = sqrt(G*Geff)/2;
    OmeC=sqrt(c3/c4)*OmeT;
    aa=((3*Gc+nh*rh)*G+nh*rh*(Gc + nh*rh))/(2*(G+1)*(Gc + nh*rh));
    % aa=((3*Gc+Wh)/2*G+Wh*(Gc+Wh)/2)/(G+Gc+Wh);
    OmeT1 =sqrt(aa*G*Geff)/2;
    P43(1,i)=(OmeT-OmeC)/(OmeC);
    Pstar(1,i)=(OmeT1-OmeC)/(OmeC);
end
end
