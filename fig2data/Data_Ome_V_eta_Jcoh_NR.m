% parallel_scan.m
% 并行参数扫描（RWA 与 Non-RWA），可直接运行
clear; close all; clc
tic

%% ---------- basis & operators ----------
e1=[1;0;0;0]; e2=[0;1;0;0]; e3=[0;0;1;0]; e4=[0;0;0;1];
P1=e1*e1'; P2=e2*e2'; P3=e3*e3'; P4=e4*e4';
S14=e1*e4'; S14d=S14';
S34=e3*e4'; S34d=S34';
S12=e1*e2'; S12d=S12';
S23=e2*e3'; S23d=S23';

%% ---------- static parameters ----------
E1 = 0; E2 = 1451; E3 = 13229; E4 = 14840;
H0 = diag([E1,E2,E3,E4]);

nh=10000; rh = 0.0016; rc=18; Gc=200;
nc=0; Nc=0;
phi=0; Delta=0; omega43=(E4-E3);

p.P1=P1; p.P2=P2; p.P3=P3; p.P4=P4;
p.S14=S14; p.S14d=S14d; p.S34=S34; p.S34d=S34d;
p.S12=S12; p.S12d=S12d; p.S23=S23; p.S23d=S23d;
p.rh=rh; p.rc=rc; p.Gc=Gc; p.nh=nh; p.nc=nc; p.Nc=Nc;
p.phi=phi; p.Delta=Delta; p.H0=H0; p.omega43=omega43;

% initial state
rho0 = P1; y0 = rho0(:);

%% ---------- scan ranges ----------
Omelist = linspace(0, 1, 100)*1611;     % 列方向（Rabi 振幅）
Glist   = logspace(log10(1e1), log10(1e6), 100); % 行方向（Γ）

wdrive = p.omega43 + p.Delta;
T = 2*pi/abs(wdrive);

nPeriods_relax = 300;                 % 更长松弛
nPeriods_avg   = 20;                  % 稳态后多周期平均（用于设置积分跨度）
PtsPerPeriod   = 400;                 % 每周期采样密度

Tend_all  = (nPeriods_relax+nPeriods_avg)*T;
Nt_all    = (nPeriods_relax+nPeriods_avg)*PtsPerPeriod;
tspan_all = linspace(0, Tend_all, Nt_all); 

% 为 RWA/Non-RWA 统一采用松弛阶段跨度（与原脚本一致）
Tend_relax = nPeriods_relax*T;
tspan_relax = linspace(0, Tend_relax, nPeriods_relax*40);
opts  = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',T/PtsPerPeriod);

nO = numel(Omelist);
nG = numel(Glist);

% 线性索引展开（parfor 友好）
[IO, JG] = ndgrid(1:nO, 1:nG);
IO = IO(:); JG = JG(:);
Ngrid = numel(IO);

% 线性结果容器（长度 = Ngrid）
RWA_eta_lin   = zeros(Ngrid,1);
NR_eta_lin    = zeros(Ngrid,1);

%% ---------- 并行池 ----------
try
    gcp('nocreate');
    if isempty(ans) %#ok<NOANS>
        parpool('threads'); % 若 R2023a+ 可选 threads 池；否则用 parpool('local')
    end
catch
    parpool('local');
end

%% ---------- 并行池(限流) ----------
% numCores   = feature('numcores');
% maxWorkers = max(1, floor(0.5 * numCores));
% 
% % 限制 MKL/BLAS 线程数，避免单个 worker 抢占全部核（如不需要可注释）
% try
%     maxNumCompThreads( max(1, ceil(maxWorkers/2)) );
% catch
%     % 某些版本可能弃用，忽略即可
% end
% 
% % 启动受限大小的本地并行池
% pool = gcp('nocreate');
% if isempty(pool) || pool.NumWorkers ~= maxWorkers
%     if ~isempty(pool), delete(pool); end
%     parpool('local', maxWorkers);
% end
%% ---------- 并行参数扫描 ----------
parfor n = 1:Ngrid
    % 局部参数
    iO = IO(n); jG = JG(n);
    p_loc = p;                         % 复制结构体，避免写时冲突
    p_loc.Omega = Omelist(iO);
    p_loc.G     = Glist(jG);

    % --------- RWA (interaction picture, autonomous) ----------
    [tRWA, yRWA] = ode15s(@(t,y) rhs_rwa(t,y,p_loc), tspan_relax, y0, opts); 
    % 取最后一个周期时间窗
    selR = tRWA >= (tspan_relax(end) - T);
    y_last_R = yRWA(selR,:);
    Nr = size(y_last_R,1);
    rho33_R = zeros(Nr,1);
    for k=1:Nr
        R = reshape(y_last_R(k,:),4,4);
        rho33_R(k) = real(R(3,3));
    end
    V = (208.5 * log(200/p_loc.G) + 11778) / 11778;
    RWA_eta_lin(n)   = p_loc.G*mean(rho33_R)*V/(0.95*p_loc.nh*p_loc.rh*(E4-E1)/11634);

    % --------- Non-RWA (lab frame, periodic) ----------
    [tNR, yNR] = ode15s(@(t,y) rhs_norwa(t,y,p_loc), tspan_relax, y0, opts);
    selN = tNR >= (tspan_relax(end) - T);
    y_last_N = yNR(selN,:);
    t_last_N = tNR(selN);
    Nn = size(y_last_N,1);
    rho33_N = zeros(Nn,1);
    for k=1:Nn
        R = reshape(y_last_N(k,:),4,4);
        rho33_N(k) = real(R(3,3));
    end
    V = (208.5 * log(200/p_loc.G) + 11778) / 11778;
    NR_eta_lin(n)   = p_loc.G*mean(rho33_N)*V/(0.95*p_loc.nh*p_loc.rh*(E4-E1)/11634);
end

%% ---------- 线性 -> 矩阵 ----------
% 注意：原图像维度为 (rows = G, cols = Omega)
RWA_eta  = reshape(RWA_eta_lin , [nO, nG]).';
NR_eta   = reshape(NR_eta_lin  , [nO, nG]).';

% 为与原脚本一致，这里转置回来：G 为行，Omega 为列
% RWA_eta  = RWA_eta.';
% NR_eta   = NR_eta.';

toc

VV = (208.0 * log(200.0 ./Glist) + 11778) /11778;

%% ---------- plot: 2x2 heatmaps ----------
figure('Color','w','Position',[100 100 1100 500]);

subplot(1,2,1);
contourf(Omelist/1611, VV, RWA_eta, 500, 'LineColor','none');
colorbar
xlim([0 4]); ylim([0.85 1]);
title('RWA \eta')

subplot(1,2,2);
contourf(Omelist/1611, VV, NR_eta, 500, 'LineColor','none');
colorbar
xlim([0 4]); ylim([0.85 1]);
title('RWA J_{coh}')


%% ---------- save ----------
% save('RWA_eta.mat', 'RWA_eta');
% save('NR_eta.mat', 'NR_eta');

%% ---------- local functions (for workers) ----------
function dy = rhs_rwa(t,y,p)
    rho = reshape(y,4,4);
    wdrive = p.omega43 + p.Delta;
    Hdrive = p.Omega*( exp(-1i*(wdrive*t + p.phi))*p.S34+ exp( 1i*(wdrive*t + p.phi))*p.S34d);
    Hlab = p.H0 + Hdrive;
    drho = -1i*(Hlab*rho - rho*Hlab);
    drho = drho + dissipators(rho,p);
    dy = drho(:);
end


function dy = rhs_norwa(t,y,p)
    rho = reshape(y,4,4);
    wdrive = p.omega43 + p.Delta;
    Hdrive = 2*p.Omega*cos(wdrive*t + p.phi) * (p.S34 + p.S34d);
    Hlab   = p.H0 + Hdrive;
    drho = -1i*(Hlab*rho - rho*Hlab);
    drho = drho + dissipators(rho,p);
    dy = drho(:);
end


function L = D(rho,A)
    L = A*rho*A' - 0.5*(A'*A*rho + rho*A'*A);
end

function drho = dissipators(rho,p)
    drho = zeros(4,4);
    drho = drho + (p.rh)*((p.nh+1)*D(rho,p.S14) + p.nh*D(rho,p.S14d));
    drho = drho + (p.rc)*((p.nc+1)*D(rho,p.S34) + p.nc*D(rho,p.S34d));
    drho = drho + (p.Gc)*((p.Nc+1)*D(rho,p.S12) + p.Nc*D(rho,p.S12d));
    drho = drho + (p.G/2)*(2*p.S23*rho*p.S23d - p.S23d*p.S23*rho - rho*p.S23d*p.S23);
end
