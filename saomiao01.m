clear all;
clc;
addpath(genpath('D:\Matlab2021b\DrosteEffect-BrewerMap-3.2.8.0'));
assert(~isempty(which('brewermap')), 'brewermap 未在路径上');

cmap = brewermap(256,'-RdYlBu');   % 需要反向就用 '-RdYlBu'
set(groot,'DefaultFigureColormap',cmap);  % 之后新建的所有 figure 默认用此色表
tic
Glist = logspace(log10(1e1), log10(1e7), 200);
Omelist = linspace(0, 0.5, 200)*1611;
% rh = 0.05;
% rc = 140;
% Gc= 200;
% nh=60000;

rh = 0.0016;
rc = 18;
Gc= 200;
nh=10000;


nc = 0;
a3 = 0.5 * (Glist + rc*nc);
a4 = 0.5 * (rh*(nh+1) + rc*(nc+1));
Ocrit =sqrt(a3.*a4);
% Ocrit = sqrt(Glist*(rc+Gc));
Xcrit = Ocrit / 1611;
Ycrit = (208.5*log(Gc./Glist) + 11778.0) / 11778;
[~,ord] = sort(Xcrit);
Xcrit = Xcrit(ord);  Ycrit = Ycrit(ord);
EPline = abs(88001 -200*Glist) / 800/11778;
V = (208.5*log(Gc./Glist) + 11778.0)/11778;
[eta,P43] = fun01(Glist, Omelist,rh,rc,Gc,nh);
[X, Y] = meshgrid(Omelist/1611, V); % x, y 是原始坐标
toc
%%
close all
h = figure('Position', [-1600,150,1000,500]);
ax1 = axes('Position', [0.1,0.2,0.35,0.67]); 
contourf(X, Y, eta, 200, 'LineColor', 'none');hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2);
cb1 = colorbar(ax1,'Location','northoutside');
xlim([min(Omelist/1611) max(Omelist/1611)])
ylim([0.85 1.05])

ax2 = axes('Position', [0.6,0.2,0.35,0.67]); 
contourf(X, Y, P43, 200, 'LineColor', 'none');hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2);
cb2 = colorbar(ax2,'Location','northoutside');
xlim([min(Omelist/1611) max(Omelist/1611)])
ylim([min(V) max(V)])


function [eta,P43] = fun01(Glist, Omelist,rh,rc,Gc,nh)
Eac=1611;
Eab=11778;
n=length(Glist);
m=length(Omelist);

eta=zeros(m,n);
jj=zeros(m,n);
PP=zeros(m,n);
Pin=zeros(m,n);
Psun=zeros(m,n);
P43=zeros(m,n);
for i =1:n
    Omeac=Omelist(i);
    for j = 1:m
        G=Glist(j);
        V = (208.5*log(Gc/G) + 11778.0)/11778;
        jj(j,i)=(1200000.0*G*(40.0*Omeac^2 + 1400.0*G + 616007.0))/(38000200.0*G^2 + 400000.0*G*Omeac^2 + 18400278001.0*G + 144000800.0*Omeac^2 + 739208400000.0);
        eta(j,i)=(7.09e-5*G*Gc*(208.0*log(Gc/G) + 1.18e+4)*(G*rc + rc*rh + 4.0*Omeac^2 + rc^2 + nh*rc*rh))/(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2);
        P43(j,i)=(G*Gc*(208.5*log(Gc/G) + 11778.0)*(4.0*Omeac^2*nh*rh + nh^2*rc*rh^2 + nh*rc^2*rh + nh*rc*rh^2 + G*nh*rc*rh))/(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2);
    end
end
end
