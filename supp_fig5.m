% parallel_scan.m
clear all;close all;clc;
nh=10000;
rh = 0.0018;
rc = 18;
Gc= 200;
E41=14856;
etaC=0.95;
tic
V1=0.92;
V2=0.94;
V3=0.96;
V4=0.98;
Omelist  = linspace(0, 1,100)*1611;
G1 = 200*exp(11778*(1-V1)/208.5);
G2 = 200*exp(11778*(1-V2)/208.5);
G3 = 200*exp(11778*(1-V3)/208.5);
G4 = 200*exp(11778*(1-V4)/208.5);
OmeT1 = sqrt(0.5*G1*(0.5*(rh*(nh+1) + rc)))/1611;
OmeT2 = sqrt(0.5*G2*(0.5*(rh*(nh+1) + rc)))/1611;
OmeT3 = sqrt(0.5*G3*(0.5*(rh*(nh+1) + rc)))/1611;
OmeT4 = sqrt(0.5*G4*(0.5*(rh*(nh+1) + rc)))/1611;
[etalist1,detadO2list1,rho43list1,Jcohlist1] = fun01(V1,nh,rh,rc,Gc,Omelist);
[etalist2,detadO2list2,rho43list2,Jcohlist2] = fun01(V2,nh,rh,rc,Gc,Omelist);
[etalist3,detadO2list3,rho43list3,Jcohlist3] = fun01(V3,nh,rh,rc,Gc,Omelist);
[etalist4,detadO2list4,rho43list4,Jcohlist4] = fun01(V4,nh,rh,rc,Gc,Omelist);
Omelist  = linspace(0, 1,100);
etamin1 = G1*Gc*11778*V1*rc/(E41*etaC*(G1*Gc*(rc+2*nh*rh)+(G1+Gc)*nh*rc*rh));
etamax1 = G1*Gc*11778*V1/(E41*etaC*(G1*(Gc+nh*rh)+3*nh*rh*Gc));

etamin2 = G2*Gc*11778*V2*rc/(E41*etaC*(G2*Gc*(rc+2*nh*rh)+(G2+Gc)*nh*rc*rh));
etamax2 = G2*Gc*11778*V2/(E41*etaC*(G2*(Gc+nh*rh)+3*nh*rh*Gc));

etamin3 = G3*Gc*11778*V3*rc/(E41*etaC*(G3*Gc*(rc+2*nh*rh)+(G3+Gc)*nh*rc*rh));
etamax3 = G3*Gc*11778*V3/(E41*etaC*(G3*(Gc+nh*rh)+3*nh*rh*Gc));

etamin4 = G4*Gc*11778*V4*rc/(E41*etaC*(G4*Gc*(rc+2*nh*rh)+(G4+Gc)*nh*rc*rh));
etamax4 = G4*Gc*11778*V4/(E41*etaC*(G4*(Gc+nh*rh)+3*nh*rh*Gc));
toc
%% 
clc
close all
brcolor=[237,249,249]/255;
excolor=[254,245,232]/255;
FaceA=1;
textbrcolor=[0 0 0];
textexcolor=[0 0 0];
lft=0.1;
hgap = 0.08;  % 水平间隙
width = 0.15;  % 子图宽度
height = 0.33; % 子图高度
bot_bot=0.15;
top_bot=0.54;

h = figure('Position', [-1900,150,2000,800]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = axes('Position', [lft, top_bot, width, height]); 
hold(ax1,'on');
yyaxis left
patch(ax1, [0 OmeT1 OmeT1 0], [8 8 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax1, [OmeT1 1 1 OmeT1], [8 8 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax1, OmeT1*ones(1,20), linspace(5,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax1, linspace(0.001,50,20),100*etamin1*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax1, linspace(0.001,50,20),100*etamax1*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h1=plot(Omelist, 100*etalist1, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax1 = gca;
ax1.YColor = [0 0 0];
ylim(ax1, [20 80]);
yticks(ax1, [20 100*etamin1 50 100*etamax1 80]);
yticklabels(ax1, {'20%','','50%','','80%'});
text(ax1, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=0.92$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)

yyaxis right
h2=plot(Omelist, 1e6*detadO2list1, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax1, [-2 6]);
yticks(ax1, [-2 0 6]);
yticklabels(ax1, {'-2','0','6'});
text(1.085, 0.55,'($10^{-6}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax1 = gca;
ax1.YColor = [0 0 0];
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 18);
xlim(ax1, [0 1]);
xticks([0 OmeT1 0.5 1]);
xticklabels({''});
text(-0.2, 0.13, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, 0.88, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
lgd1=legend([h1,h2], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd1.Position = [0.16, 0.65, 0.1, 0.1];  
lgd1.Box = 'off';
lgd1.FontSize = 22;
text(-0.2,1.2, '(a)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax2 = axes('Position', [lft+width+hgap, top_bot, width, height]); 
hold(ax2,'on');
yyaxis left
patch(ax2, [0 OmeT2 OmeT2 0], [8 8 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax2, [OmeT2 1 1 OmeT2], [8 8 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax2, OmeT2*ones(1,20), linspace(5,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax2, linspace(0.001,50,20),100*etamin2*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax2, linspace(0.001,50,20),100*etamax2*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h3=plot(Omelist, 100*etalist2, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax2 = gca;
ax2.YColor = [0 0 0];
ylim(ax2, [20 80]);
yticks(ax2, [20 100*etamin2 50 100*etamax2 80]);
yticklabels(ax2, {'20%','','50%','','80%'});
text(ax2, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=0.94$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)

yyaxis right
h4=plot(Omelist, 1e6*detadO2list2, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax2, [-5 15]);
yticks(ax2, [-5 0 15]);
yticklabels(ax2, {'-5','0','15'});
text(1.085, 0.55,'($10^{-6}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax2 = gca;
ax2.YColor = [0 0 0];
ax2.Box      = 'on';
ax2.Layer    = 'top';
ax2.LineWidth = 2;
ax2.TickDir  = 'in';
set(ax2, 'FontSize', 18);
xlim(ax2, [0 1]);
xticks([0 OmeT2 0.5 1.0]);
xticklabels({''});
text(-0.2, 0.14, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, 0.87, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
lgd2=legend([h3,h4], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd2.Position = [0.39, 0.65, 0.1, 0.1];  
lgd2.Box = 'off';
lgd2.FontSize = 22;
text(-0.2,1.2, '(b)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax3 = axes('Position', [lft+2*width+2*hgap, top_bot, width, height]); 
hold(ax3,'on');
yyaxis left
patch(ax3, [0 OmeT3 OmeT3 0], [8 8 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax3, [OmeT3 1 1 OmeT3], [8 8 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax3, OmeT3*ones(1,20), linspace(5,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax3, linspace(0.001,50,20),100*etamin3*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax3, linspace(0.001,50,20),100*etamax3*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h5=plot(Omelist, 100*etalist3, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax3 = gca;
ax3.YColor = [0 0 0];
ylim(ax3, [20 80]);
yticks(ax3, [20 100*etamin3 50 100*etamax3 80]);
yticklabels(ax3, {'20%','','50%','','80%'});
text(ax3, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=0.96$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)


yyaxis right
h6=plot(Omelist, 1e5*detadO2list3, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax3, [-2 4]);
yticks(ax3, [-2 0 4]);
yticklabels(ax3, {'-2','0','4'});
text(1.085, 0.55,'($10^{-5}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax3 = gca;
ax3.YColor = [0 0 0];
ax3.Box      = 'on';
ax3.Layer    = 'top';
ax3.LineWidth = 2;
ax3.TickDir  = 'in';
set(ax3, 'FontSize', 18);
xlim(ax3, [0 1]);
xticks([0 OmeT3 0.5 1.0]);
xticklabels({''});
text(-0.2, 0.15, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, 0.88, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2,1.2, '(c)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);
lgd3=legend([h5,h6], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd3.Position = [0.62, 0.65, 0.1, 0.1];  
lgd3.Box = 'off';
lgd3.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax4 = axes('Position', [lft+3*width+3*hgap, top_bot, width, height]); 
hold(ax4,'on');
yyaxis left
patch(ax4, [0 OmeT4 OmeT4 0], [8 8 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax4, [OmeT4 1 1 OmeT4], [8 8 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax4, OmeT4*ones(1,20), linspace(5,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax4, linspace(0.001,50,20),100*etamin4*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax4, linspace(0.001,50,20),102*etamax4*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h5=plot(Omelist, 100*etalist4, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax4 = gca;
ax4.YColor = [0 0 0];
ylim(ax4, [20 80]);
yticks(ax4, [20 100*etamin4 50 102*etamax4 80]);
yticklabels(ax4, {'20%','','50%','','80%'});
text(ax4, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=0.98$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)


yyaxis right
h6=plot(Omelist, 1e5*detadO2list4, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax4, [-5 15]);
yticks(ax4, [-5 0 15]);
yticklabels(ax4, {'-5','0','15'});
text(1.085, 0.55,'($10^{-5}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax4 = gca;
ax4.YColor = [0 0 0];
ax4.Box      = 'on';
ax4.Layer    = 'top';
ax4.LineWidth = 2;
ax4.TickDir  = 'in';
set(ax4, 'FontSize', 18);
xlim(ax4, [0 1]);
xticks([0 OmeT4 0.5 1.0]);
xticklabels({''});
text(-0.2, 0.17, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, 0.7, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2,1.2, '(d)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);
lgd4=legend([h5,h6], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd4.Position = [0.85, 0.65, 0.1, 0.1];  
lgd4.Box = 'off';
lgd4.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax5 = axes('Position', [lft, bot_bot, width, height]); 
hold(ax5,'on');
yyaxis left
patch(ax5, [0 OmeT1 OmeT1 0], [-1 -1 0.12 0.12],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax5, [OmeT1 4 4 OmeT1], [-1 -1 0.12 0.12],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax5, OmeT1*ones(1,20), linspace(-0.1,0.1,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h7=plot(Omelist, rho43list1, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax5.YColor = [0 0 0];
ylim(ax5, [-0.0005 0.01]);
yticks(ax5, [0 0.005 0.01]);
yticklabels(ax5, {'0', '5','10'});
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax5, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(ax5, 0.03, 1.04, '$\times 10^{-3}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h8=plot(Omelist, Jcohlist1, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax5, [-0.8 16]);
yticks(ax5, [0 16]);
yticklabels(ax5, {'0','16'});
ax.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax5 = gca;
ax5.YColor = [0 0 0];
ax5.Box      = 'on';
ax5.Layer    = 'top';
ax5.LineWidth = 2;
ax5.TickDir  = 'in';
set(ax5, 'FontSize', 18);
xlim(ax5, [0 1]);
xticks([0 OmeT1 0.5 1.0]);
xticklabels({'0','','0.5','1.0'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(0.22, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd5=legend([h7, h8], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd5.Position = [0.16, 0.2, 0.1, 0.1];  
lgd5.Box = 'off';
lgd5.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax6 = axes('Position', [lft+width+hgap, bot_bot, width, height]); 
hold(ax6,'on');
yyaxis left
patch(ax6, [0 OmeT2 OmeT2 0], [-1 -1 0.12 0.12],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax6, [OmeT2 4 4 OmeT2], [-1 -1 0.12 0.12],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax6, OmeT2*ones(1,20), linspace(0,0.12,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h9=plot(Omelist, rho43list2, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax6.YColor = [0 0 0];
ylim(ax6, [-0.001 0.02]);
yticks(ax6, [0 0.01 0.02]);
yticklabels(ax6, {'0', '1','2'});
text(ax6, 0.03, 1.04, '$\times 10^{-2}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax6, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h10=plot(Omelist, Jcohlist2, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax6, [-0.8 16]);
yticks(ax6, [-0.8 16]);
yticklabels(ax6, {'0','16'});
ax.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax6 = gca;
ax6.YColor = [0 0 0];
ax6.Box      = 'on';
ax6.Layer    = 'top';
ax6.LineWidth = 2;
ax6.TickDir  = 'in';
set(ax6, 'FontSize', 18);
xlim(ax6, [0 1]);
xticks([0 OmeT2 0.5 1.0]);
xticklabels({'0','','0.5','1.0'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd5=legend([h9, h10], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd5.Position = [0.34, 0.18, 0.1, 0.1];  
lgd5.Box = 'off';
lgd5.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax7 = axes('Position', [lft+2*width+2*hgap, bot_bot, width, height]); 
hold(ax7,'on');
yyaxis left
patch(ax7, [0 OmeT3 OmeT3 0], [-1 -1 0.12 0.12],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax7, [OmeT3 4 4 OmeT3], [-1 -1 0.12 0.12],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax7, OmeT3*ones(1,20), linspace(-0.1,0.12,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h11=plot(Omelist, rho43list3, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax7.YColor = [0 0 0];
ylim(ax7, [-0.006 0.04]);
yticks(ax7, [0 0.02 0.04]);
yticklabels(ax7, {'0', '2','4'});
text(ax7, 0.03, 1.04, '$\times 10^{-2}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax7, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h12=plot(Omelist, Jcohlist3, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax7, [-0.8 16]);
yticks(ax7, [-0.8 16]);
yticklabels(ax7, {'0','16'});
ax.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax7 = gca;
ax7.YColor = [0 0 0];
ax7.Box      = 'on';
ax7.Layer    = 'top';
ax7.LineWidth = 2;
ax7.TickDir  = 'in';
set(ax7, 'FontSize', 18);
xlim(ax7, [0 1]);
xticks([0 OmeT3 0.5 1.0]);
xticklabels({'0','','0.5','1.0'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(0.08, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd4=legend([h11, h12], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd4.Position = [0.62, 0.28, 0.1, 0.1];  
lgd4.Box = 'off';
lgd4.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
ax8 = axes('Position', [lft+3*width+3*hgap, bot_bot, width, height]); 
hold(ax8,'on');
yyaxis left
patch(ax8, [0 OmeT4 OmeT4 0], [-2 -2 0.2 0.2],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax8, [OmeT4 4 4 OmeT4], [-2 -2 0.2 0.2],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax8, OmeT4*ones(1,20), linspace(0,0.12,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h11=plot(Omelist, rho43list4, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax8.YColor = [0 0 0];
ylim(ax8, [-0.003 0.06]);
yticks(ax8, [0 0.03 0.06]);
yticklabels(ax8, {'0', '3','6'});
text(ax8, 0.03, 1.04, '$\times 10^{-2}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax8, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h12=plot(Omelist, Jcohlist4, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax8, [-0.8 16]);
yticks(ax8, [-0.8 16]);
yticklabels(ax8, {'0','16'});
ax.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax8 = gca;
ax8.YColor = [0 0 0];
ax8.Box      = 'on';
ax8.Layer    = 'top';
ax8.LineWidth = 2;
ax8.TickDir  = 'in';
set(ax8, 'FontSize', 18);

xlim(ax8, [0 1]);
xticks([0 OmeT4 0.5 1.0]);
xticklabels({'0','','0.5','1.0'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(0.04, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd4=legend([h11, h12], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd4.Position = [0.85, 0.28, 0.1, 0.1];  
lgd4.Box = 'off';
lgd4.FontSize = 22;



%% 
function [etalist,detadO2list,rho43_1D,Jcoh_1D] = fun01(V,nh,rh,rc,Gc,Omelist)
G= 200*exp(11778*(1-V)/208.5);
nO = numel(Omelist);
etalist   = zeros(1,nO);
detadO2list  = zeros(1,nO);
for iO = 1:nO
    Omeac = Omelist(iO);
    etalist(iO) = (0.088569566079981860952866819714877*G*((417*log(200/G))/2 + 11778)*(1250.0*Omeac^2 + 5625.0*G + 191259.0))/(20094375.0*G^2 + 1687500.0*G*Omeac^2 + 795740901.0*G + 75002500.0*Omeac^2 + 3825180000.0);
    detadO2list(iO)=(7812500*G*((417*log(200/G))/2 + 11778))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)) - (15625000*G*Omeac*(150005000*Omeac + 3375000*G*Omeac)*((417*log(200/G))/2 + 11778))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)^2) + (6250*G*(150005000*Omeac + 3375000*G*Omeac)^2*((417*log(200/G))/2 + 11778)*(1250*Omeac^2 + 5625*G + 191259))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)^3) - (3125*G*(3375000*G + 150005000)*((417*log(200/G))/2 + 11778)*(1250*Omeac^2 + 5625*G + 191259))/(35283*(20094375*G^2 + 1687500*G*Omeac^2 + 795740901*G + 75002500*Omeac^2 + 3825180000)^2);
    rho43=(1.25e+7*(18.0*Omeac - 1.0*G*Omeac))/(G^2*2.01e+7i + G*Omeac^2*1.69e+6i + G*7.96e+8i + Omeac^2*7.5e+7i + 3.83e+9i);
    rho43_1D(iO)=abs(rho43);
    Jcoh_1D(iO) =2*Omeac*abs(rho43);
end
end