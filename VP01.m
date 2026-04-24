clear; clc; close all 
syms rho44 rho43 rho42 rho41 rho34 rho33 rho32 rho31 ...
     rho24 rho23 rho22 rho21 rho14 rho13 rho12 rho11
syms G Omeac Omeab rh rc Gc kB Ts Ta nh nc Nc E32 E41 Del del V real
nh = sym(10000);
rh = sym(0.0016);
rc = sym(18);
Gc = sym(200);
kB = sym(0.6950);
Ts = sym(6000);
Ta = sym(300);
E32 = sym(11778);
E41 = sym(14856);
nc  = sym(0);
Nc  = sym(0);
del=0;
Del=0;
V=0.95;
G = 200*exp(11778*(1-V)/208.5);
% G=75.5;
OmeT = sqrt(0.5 * (G + rc*nc)*0.5 * (rh*(nh+1) + rc*(nc+1)));

c3=(G+rc+nh*rh)*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh);
c4=G*(rc+rh*nh)*(G*Gc+G*nh*rh+3*Gc*nh*rh);
OmeC=sqrt(c3/c4)*OmeT;
% G=3.3704e+03;
% Omeac= 0;
% Omeac= OmeC;
% Omeac= OmeT;
% Omeac= 1E20*OmeC;
Omeab=0;
etaC=0.95;
eqns = [
Gc*rho22 + rh*rho44*(nh + 1) - nh*rh*rho11==0
- (Gc*rho12)/2 - (nh*rh*rho12)/2==0
Omeac*rho14*1i - (G*rho13)/2 - (nh*rh*rho13)/2==0
Omeac*rho13*1i - (rc*rho14)/2 - (rh*rho14*(nh + 1))/2 - (nh*rh*rho14)/2==0
- (Gc*rho21)/2 - (nh*rh*rho21)/2==0
G*rho33 - Gc*rho22==0
Omeac*rho24*1i - (Gc*rho23)/2 - (G*rho23)/2==0
Omeac*rho23*1i - (Gc*rho24)/2 - (rc*rho24)/2 - (rh*rho24*(nh + 1))/2==0
- (G*rho31)/2 - Omeac*rho41*1i - (nh*rh*rho31)/2==0
- (G*rho32)/2 - (Gc*rho32)/2 - Omeac*rho42*1i==0
Omeac*rho34*1i - G*rho33 - Omeac*rho43*1i + rc*rho44==0
Omeac*rho33*1i - (G*rho34)/2 - Omeac*rho44*1i - (rc*rho34)/2 - (rh*rho34*(nh + 1))/2==0
- Omeac*rho31*1i - (rc*rho41)/2 - (rh*rho41*(nh + 1))/2 - (nh*rh*rho41)/2==0
- (Gc*rho42)/2 - Omeac*rho32*1i - (rc*rho42)/2 - (rh*rho42*(nh + 1))/2==0
Omeac*rho44*1i - Omeac*rho33*1i - (G*rho43)/2 - (rc*rho43)/2 - (rh*rho43*(nh + 1))/2==0
Omeac*rho43*1i - Omeac*rho34*1i - rc*rho44 - rh*rho44*(nh + 1) + nh*rh*rho11==0
rho44 + rho33 + rho22 + rho11 == 1
];

% 用 equationsToMatrix 抽取线性系统（避免 solve 的变量检查路径）
vars = [rho44; rho43; rho42; rho41; rho34; rho33; rho32; rho31; ...
        rho24; rho23; rho22; rho21; rho14; rho13; rho12; rho11];

[A,b] = equationsToMatrix(eqns, vars);   % 这里若成功，即证明 vars 的类型合法
x = A \ b;                                % 线性精确解（symbolic）
names = {'rho44','rho43','rho42','rho41','rho34','rho33','rho32','rho31', ...
         'rho24','rho23','rho22','rho21','rho14','rho13','rho12','rho11'};
S = cell2struct(cellfun(@(k) x(k), num2cell(1:numel(names)), 'UniformOutput', false), names, 2);
if isstruct(S)
    rho = S;
else
    % 若解为向量 x，则需按与 vars 一致的顺序映射
    names = {'rho44','rho43','rho42','rho41','rho34','rho33','rho32','rho31',...
             'rho24','rho23','rho22','rho21','rho14','rho13','rho12','rho11'};
    rho = cell2struct(num2cell(S), names, 2);
end
rhoaa_value = vpa(rho.rho44,3);
rhoac_value = vpa(rho.rho43,3);
rhoav_value = vpa(rho.rho42,3);
rhoab_value = vpa(rho.rho41,3);
rhoca_value = vpa(rho.rho34,3);
rhocc_value = vpa(rho.rho33,3);
rhocv_value = vpa(rho.rho32,3);
rhocb_value = vpa(rho.rho31,3);
rhova_value = vpa(rho.rho24,3);
rhovc_value = vpa(rho.rho23,3);
rhovv_value = vpa(rho.rho22,3);
rhovb_value = vpa(rho.rho21,3);
rhoba_value = vpa(rho.rho14,3);
rhobc_value = vpa(rho.rho13,3);
rhobv_value = vpa(rho.rho12,3);
rhobb_value = vpa(rho.rho11,3);

%% 计算并输出物理量
clc
VV = E32 + kB*Ta*log(rho.rho33/rho.rho22);
jj = vpa(G*rho.rho33,3);

% disp(vpa(jj*VV/(etaC*E41*nh*rh),3))

% disp(rhoaa_value)
% disp(rhocc_value)
% % disp(rhoac_value)
% disp(simplify(jj*VV/(etaC*E41*nh*rh)))
disp(simplify(VV))
% disp(simplify(jj*VV))

% disp(rhocc_value*G/Gc)
% disp(rhobb_value)

% disp(rhoac_value)
% disp(rhoca_value)
% disp(rhoav_value)
% disp(rhoab_value)
% disp(rhocv_value)
% disp(rhocb_value)
% disp(rhovb_value)
% disp(simplify(rhoaa_value))
% disp(simplify(rhocc_value))
% disp(simplify(rhoac_value))
% disp(simplify(rhovv_value))
% disp(simplify(rhobb_value))


% disp(imag(rhoac_value))
% disp(abs(rhoac_value))
% disp(simplify(rhocc_value))
% disp(simplify(rhoaa_value-rhocc_value))
% disp(simplify(rhovv_value))
% disp(simplify(imag(rhovb_value)))
% disp(simplify(rhobb_value))
% disp(simplify(rhoaa_value+rhocc_value+rhovv_value+rhobb_value))
% disp(vpa(jj*VV))

% disp(vpa(simplify(jj*V/(etaC*E41*nh*rh)),3))
% disp(vpa(simplify(G*rhocc_value*VV/(etaC*E41*nh*rh)),3))