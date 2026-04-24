clear all;
clc;
nh=10000;
rh = 0.0018;
rc = 18;
Gc= 200;
kB=0.6950;
Ts=6000;
Ta=300;
Vlist = linspace(0.8,1.0,100)*11778;
Omelist = logspace(-2, log10(1), 200)*11778;

q1=1;
q2=1.0;
% q3=0.93;
q3=1;
q4=1.0;
OmeCTlist=zeros(1,length(Vlist));
for i=1:length(Vlist)
    V=Vlist(i);
    G = 200*exp(11778*(1-V/11778)/208.5);
    c3=(G+rc+nh*rh)*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh);
    c4=G*(rc+rh*nh)*(G*Gc+G*nh*rh+3*Gc*nh*rh);
    OmeT = sqrt(0.5*G*(0.5*(rh*(nh+1) + rc)));
    OmeC=sqrt(c3/c4)*OmeT;
    rh=q1*rh;
    % rc=q2*rc;
    % G=q3*G;
    % Gc=q4*Gc;
    Geff=(rh*(nh+1) + rc);
    % Geff1=(Geff+G)/2+0.5*sqrt((Geff-G)^2+4*Geff*G*(c3/c4-1));
    % G1=(Geff+G)/2-0.5*sqrt((Geff-G)^2+4*Geff*G*(c3/c4-1));
    % OmeT1 = q3*sqrt(0.5*G*(0.5*abs(Geff)));
    % aa = (19*G + 8360)/(22*G + 7920);
    % % Geff1=Geff*(1-(rh*(nh+1)/(rh*(nh+1)+Gc)*G/(G+Geff)));
    % % aa = (Geff*(G + Geff))/(G*(Gc + nh*rh) + 3*nh);
    % aa = ((rh*nh + rc)*(G + (rh*nh + rc)))/(G*(Gc + nh*rh) + 3*nh);
    % aa = ((rh*nh + rc)*(G + (rh*nh + rc)))/(G*(Gc + nh*rh) + 3*nh);
    aa=((3*Gc+nh*rh)*G+nh*rh*(Gc + nh*rh))/(2*(G+1)*(Gc + nh*rh));
    OmeT1 = q3*sqrt(aa*0.5*G*(0.5*(rh*(nh+1) + rc)));
    % OmeT1=0.5*sqrt(G*Geff1);
    OmeCTlist(1,i)=OmeC/OmeT1;
    OmeCTlist1(1,i)=OmeC/OmeT;
end

%%
close all
color1=[218,243,241]/255;
color2=[150,150,150]/255;
figure('Position', [-1000,200,700,450]);
ax1= subplot(1,1,1);
patch(ax1, [0.9 1.1 1.1 0.9], [-0.5 -0.5 6 6],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax1, [1.05 1.0 1.0 1.05], [-0.5 -0.5 6 6],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
h1=plot(Vlist/11778,OmeCTlist,'LineWidth',6,'color',[238,90,57]/255);hold on;
% h2=plot(Vlist/11778,OmeCTlist1,'LineWidth',4,'color',[100,136,192]/255);hold on;
plot(linspace(0.8,0.95,20), 0.998*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(0.95*ones(1,20), linspace(0.9,0.999,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(0.95, 0.998, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [238,90,57]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.00 1.05]);
xticklabels({'0.8','0.9','0.95','1.0','1.05'});
ylim([0.97 1.03])
yticks([0.97 1.0 1.03]);
yticklabels({'0.97','1.0','1.03'});
xlabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
ylabel('$\Omega_{c}/\Omega^{\ast}_{T}$','Interpreter','latex', 'FontSize', 24); 
% lgd=legend([h1, h2], {'$\Omega^{\ast}_{T}$','$\Omega_{T}$'}, 'Interpreter', 'latex');
% lgd.Position = [0.25, 0.68, 0.1, 0.1];  
% lgd.Box = 'off';
% lgd.FontSize = 22;