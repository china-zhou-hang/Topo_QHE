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
[PP,SvN,jj] = fun01(Glist, Omelist,nh,rh,rc,Gc);
[X, Y] = meshgrid(Omelist/1611, V); % x, y 是原始坐标
toc
RWA_PP=PP/11634;
RWA_SvN=SvN;
save('RWA_PP.mat', 'RWA_PP');
save('RWA_SvN.mat', 'RWA_SvN');
% save('Jcoh.mat', 'Jcoh');
%%
close all
h = figure('Position', [-1600,50,1000,600]);
ax1 = subplot(1,2,1);
contourf(X, Y, RWA_SvN,200, 'LineColor', 'none');hold on;
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
xlim([0 0.6])
ylim([0.85 1.05])
text(0.52, 1.06, '$\eta(\%)$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18);

ax2 = subplot(1,2,2);
contourf(X, Y, RWA_PP, 200, 'LineColor', 'none');hold on;
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
xlim([0 0.6])
ylim([0.85 1.05])

% ax3 = subplot(2,2,3);
% contourf(X, Y, Jcoh, 200, 'LineColor', 'none');hold on;
% plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2);
% % plot(EPline,V, 'g--', 'LineWidth', 2);hold on;
% % set(gca, 'XScale', 'log');
% xlabel('$\Omega~(units~of~\omega_{32})$', 'Interpreter', 'latex');
% ylabel('$V~(units~of~\omega_{32})$', 'Interpreter', 'latex');
% cb3 = colorbar(ax3,'Location','northoutside');
% axpos = ax3.Position;               % 主轴位置
% cb3.Position(1) = axpos(1);          % 左对齐
% cb3.Position(3) = axpos(3)/2;          % 宽度一致
% cb3.Position(2) = axpos(2) + axpos(4) + 0.02;  % 向上偏一点
% set(ax3, 'FontSize', 18);
% ax3.LineWidth = 2;
% xlim([0 1])
% ylim([0.85 1.05])


function [PP,SvN,jj] = fun01(Glist, Omelist,nh,rh,rc,Gc)
Eac=1611;
Eab=11778;
nc=0;
Nc=0;
n=length(Glist);
m=length(Omelist);
P=zeros(m,n);
SvN=zeros(m,n);
jj=zeros(m,n);
Psun=zeros(m,n);
PP=zeros(m,n);
for i =1:n
    Omeac=Omelist(i);
    for j = 1:m
        G=Glist(j);
SS=(1.0e+4*((0.5*(625.0*G - 1.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 2.13e+4)*(log(3.25e+4*G + 18.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 625.0*G^2 + 5000.0*Omeac^2 - 1.0*G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.83e+5) - 1.0*log(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) + 8.52))/(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) - (0.5*(log(3.25e+4*G - 18.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 625.0*G^2 + 5000.0*Omeac^2 + G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.83e+5) - 1.0*log(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) + 8.52)*(625.0*G + (3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 2.13e+4))/(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2))*(625.0*G^2 + 2.13e+4*G + 2500.0*Omeac^2))/(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) - ((1.0*log(1.33e+7*G^2 + 1.56e+6*G*Omeac^2 + 4.52e+8*G + 2.5e+7*Omeac^2) - log(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9))*(1.33e+7*G^2 + 1.56e+6*G*Omeac^2 + 4.52e+8*G + 2.5e+7*Omeac^2))/(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) + (2.0e+4*((0.5*(log(3.25e+4*G - 18.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 625.0*G^2 + 5000.0*Omeac^2 + G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.83e+5) - 1.0*log(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) + 8.52)*(625.0*G - 1.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 2.13e+4))/(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) - (0.5*(log(3.25e+4*G + 18.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 625.0*G^2 + 5000.0*Omeac^2 - 1.0*G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.83e+5) - 1.0*log(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) + 8.52)*(625.0*G + (3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 2.13e+4))/(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2))*(1250.0*Omeac^2 + 5620.0*G + 1.91e+5))/(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) - (100.0*G*(1250.0*Omeac^2 + 5620.0*G + 1.91e+5)*(log(1250.0*Omeac^2 + 5620.0*G + 1.91e+5) - 1.0*log(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) + log(G) + 4.61))/(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9) - (1.56e+10*Omeac^2*(1.0*log(1.63e+8*G - 9.0e+4*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.12e+6*G^2 + 2.5e+7*Omeac^2 + 5000.0*G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 1.91e+9) - log(1.63e+8*G + 9.0e+4*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.12e+6*G^2 + 2.5e+7*Omeac^2 - 5000.0*G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 1.91e+9))*(G - 18.0))/((3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2)*(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9)) - (1.25e+7*Omeac*(Omeac*log(3.25e+4*G - 18.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 625.0*G^2 + 5000.0*Omeac^2 + G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.83e+5)*1250.0i - Omeac*log(3.25e+4*G + 18.0*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 625.0*G^2 + 5000.0*Omeac^2 - 1.0*G*(3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2) + 3.83e+5)*1250.0i)*(G - 18.0)*(1.33e+7*G^2 + 1.56e+6*G*Omeac^2 + 4.52e+8*G + 2.5e+7*Omeac^2))/((3.91e+5*G^2 + 2.66e+7*G + 6.25e+6*Omeac^2 + 4.52e+8)^(1/2)*(G^2*1.33e+7i + G*Omeac^2*1.56e+6i + G*4.52e+8i + Omeac^2*2.5e+7i)*(2.01e+7*G^2 + 1.69e+6*G*Omeac^2 + 7.96e+8*G + 7.5e+7*Omeac^2 + 3.83e+9));
   SvN(j,i)=real(SS);
jj(j,i)=(20000*G*(1250*Omeac^2 + 5625*G + 191259))/((13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2)*((20000*(1250*Omeac^2 + 5625*G + 191259))/(13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2) + (10000*(625*G^2 + 21251*G + 2500*Omeac^2))/(13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2) + (100*G*(1250*Omeac^2 + 5625*G + 191259))/(13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2) + 1));
PP(j,i)=(20000*G*((417*log(200/G))/2 + 11634)*(1250*Omeac^2 + 5625*G + 191259))/((13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2)*((20000*(1250*Omeac^2 + 5625*G + 191259))/(13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2) + (10000*(625*G^2 + 21251*G + 2500*Omeac^2))/(13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2) + (100*G*(1250*Omeac^2 + 5625*G + 191259))/(13281875*G^2 + 1562500*G*Omeac^2 + 451605001*G + 25002500*Omeac^2) + 1));
    end
end
end
