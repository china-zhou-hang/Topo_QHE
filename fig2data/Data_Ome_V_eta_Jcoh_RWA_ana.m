clear all;
clc;
addpath(genpath('D:\Matlab2021b\DrosteEffect-BrewerMap-3.2.8.0'));
assert(~isempty(which('brewermap')), 'brewermap 未在路径上');

cmap = brewermap(256,'-RdYlBu');   % 需要反向就用 '-RdYlBu'
set(groot,'DefaultFigureColormap',cmap);  % 之后新建的所有 figure 默认用此色表
tic
Glist = logspace(log10(1e1), log10(1e6), 100);
Omelist = linspace(0, 0.6, 100)*1611;
% Omelist = logspace(-2, log10(30), 100)*1611;
rh = 0.0018;
rc = 18;
Gc= 200;
nh=10000;
nc = 0;
a3 = 0.5 * (Glist + rc*nc);
a4 = 0.5 * (rh*(nh+1) + rc*(nc+1));
Ocrit = sqrt(a3 .* a4);
% Ocrit = sqrt(Glist*(-rc+Gc));
Xcrit = Ocrit / 1611;
Ycrit = (208.5*log(Gc./Glist) + 11778.0) / 11778;
[~,ord] = sort(Xcrit);
Xcrit = Xcrit(ord);  Ycrit = Ycrit(ord);
EPline = abs(88001 -200*Glist) / 800/11778;
V = (208.5*log(Gc./Glist) + 11778.0)/11778;
[eta,RWA_rho43,Jcoh] = fun01(Glist, Omelist,nh,rh,rc,Gc);
[X, Y] = meshgrid(Omelist/1611, V); % x, y 是原始坐标
toc
RWA_eta=eta;
save('RWA_eta.mat', 'RWA_eta');
save('RWA_rho43.mat', 'RWA_rho43');
% save('Jcoh.mat', 'Jcoh');
%%
close all
h = figure('Position', [-1600,50,1000,1000]);
ax1 = subplot(2,2,1);
contourf(X, Y, eta,200, 'LineColor', 'none');hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2);
% plot(EPline,V, 'g--', 'LineWidth', 2);hold on;
% set(gca, 'XScale', 'log');
xlabel('$\Omega~(units~of~\omega_{32})$', 'Interpreter', 'latex');
ylabel('$V~(units~of~\omega_{32})$', 'Interpreter', 'latex');
cb1 = colorbar(ax1,'Location','northoutside');
axpos = ax1.Position;               % 主轴位置
cb1.Position(1) = axpos(1);          % 左对齐
cb1.Position(3) = axpos(3)/2;          % 宽度一致
cb1.Position(2) = axpos(2) + axpos(4) + 0.02;  % 向上偏一点
set(ax1, 'FontSize', 18);
ax1.LineWidth = 2;
xlim([0 1])
ylim([0.85 1.05])
text(0.52, 1.06, '$\eta(\%)$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18);

ax2 = subplot(2,2,2);
contourf(X, Y, 2*RWA_rho43, 200, 'LineColor', 'none');hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2);
% plot(EPline,V, 'g--', 'LineWidth', 2);hold on;
% set(gca, 'XScale', 'log');
xlabel('$\Omega~(units~of~\omega_{32})$', 'Interpreter', 'latex');
ylabel('$V~(units~of~\omega_{32})$', 'Interpreter', 'latex');
cb2 = colorbar(ax2,'Location','northoutside');
axpos = ax2.Position;               % 主轴位置
cb2.Position(1) = axpos(1);          % 左对齐
cb2.Position(3) = axpos(3)/2;          % 宽度一致
cb2.Position(2) = axpos(2) + axpos(4) + 0.02;  % 向上偏一点
set(ax2, 'FontSize', 18);
ax2.LineWidth = 2;
xlim([0 1])
ylim([0.85 1.05])

ax3 = subplot(2,2,3);
contourf(X, Y, Jcoh, 200, 'LineColor', 'none');hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2);
% plot(EPline,V, 'g--', 'LineWidth', 2);hold on;
% set(gca, 'XScale', 'log');
xlabel('$\Omega~(units~of~\omega_{32})$', 'Interpreter', 'latex');
ylabel('$V~(units~of~\omega_{32})$', 'Interpreter', 'latex');
cb3 = colorbar(ax3,'Location','northoutside');
axpos = ax3.Position;               % 主轴位置
cb3.Position(1) = axpos(1);          % 左对齐
cb3.Position(3) = axpos(3)/2;          % 宽度一致
cb3.Position(2) = axpos(2) + axpos(4) + 0.02;  % 向上偏一点
set(ax3, 'FontSize', 18);
ax3.LineWidth = 2;
xlim([0 1])
ylim([0.85 1.05])


function [eta,RWA_rho43,Jcoh] = fun01(Glist, Omelist,nh,rh,rc,Gc)
Eac=1611;
Eab=11778;
nc=0;
Nc=0;
n=length(Glist);
m=length(Omelist);
reta=zeros(m,n);
eta=zeros(m,n);
jj=zeros(m,n);
RWA_rho43=zeros(m,n);
Jcoh=zeros(m,n);
for i =1:n
    Omeac=Omelist(i);
    for j = 1:m
        G=Glist(j);
        V = (208.5*log(Gc/G) + 11778.0);
        p44=(Gc*G^2*nh*rh + Gc*G*nh^2*rh^2 + Gc*G*nh*rh^2 + Gc*rc*G*nh*rh + 4.0*Gc*Omeac^2*nh*rh)/(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2);
        p33=(Gc*(4.0*Omeac^2*nh*rh + nh^2*rc*rh^2 + nh*rc^2*rh + nh*rc*rh^2 + G*nh*rc*rh))/(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2);
        r43=-(2.0*(G*Gc*Omeac*nh*rh - 1.0*Gc*Omeac*nh*rc*rh))/(G^2*nh*rc*rh*1.0i + Gc*G^2*nh*rh*2.0i + Gc*G^2*rc*1.0i + Gc*G^2*rh*1.0i + G*Omeac^2*nh*rh*4.0i + Gc*G*Omeac^2*4.0i + G*nh^2*rc*rh^2*1.0i + Gc*G*nh^2*rh^2*2.0i + G*nh*rc^2*rh*1.0i + G*nh*rc*rh^2*1.0i + Gc*G*nh*rc*rh*4.0i + Gc*G*nh*rh^2*3.0i + Gc*G*rc^2*1.0i + Gc*G*rc*rh*2.0i + Gc*G*rh^2*1.0i + Gc*Omeac^2*nh*rh*12.0i + Gc*Omeac^2*rh*4.0i + Gc*nh^2*rc*rh^2*1.0i + Gc*nh*rc^2*rh*1.0i + Gc*nh*rc*rh^2*1.0i);
        eta(j,i)=(5*G*Gc*((417*log(Gc/G))/2 + 11778)*(G*rc + rc*rh + 4.0*Omeac^2 + rc^2 + nh*rc*rh))/(70566*(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2));
        RWA_rho43(j,i)=abs(r43);
        Jcoh(j,i)=2*Omeac*imag(r43);
    end
end
end
