%%
%Andreas Stadelmann 01525173
%Cristian Avram     01304470

clear all   % Es werden alle Variablen gelöscht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zurückgesetzt

%Eingangsgroeßen
syms uGSM Mext

%Zustandsgroeßen
syms iGSM phiGSMP wGSM wP

%Parameter
syms LGSM RGSM kGSM JGSM dcGSM dvGSM JP dcP dvP dqP cGSMP dGSMP
p_syms = {LGSM, RGSM, kGSM, JGSM, dcGSM, dvGSM, JP, dcP, dvP, dqP, cGSMP, dGSMP};

LGSM_P = 1.4e-3;
RGSM_P = 0.46;
kGSM_P = 0.1;
JGSM_P = 12.4e-3;
dcGSM_P = 0.152;
dvGSM_P = 1.8e-3;
JP_P = 32.5e-3;
dcP_P = 0.169;
dvP_P = 2.7e-3;
dqP_P = 1e-4;
cGSMP_P = 0.6822;
dGSMP_P = 1e-5;

p = [LGSM_P, RGSM_P, kGSM_P, JGSM_P, dcGSM_P, dvGSM_P, JP_P, dcP_P, dvP_P, dqP_P, cGSMP_P, dGSMP_P]

%Gleichungen 
MP = dcP + dvP*wP+dqP*wP^2+Mext;
MrGSM = dcGSM + dvGSM*wGSM;
Mkopp = (wGSM-wP)*dGSMP + phiGSMP*cGSMP;
MGSM = kGSM*iGSM;

%% 1 System reduzieren

xM = [phiGSMP; wGSM; wP];
fM = [wGSM-wP; (MGSM-MrGSM-Mkopp)/(JGSM); (Mkopp-MP)/JP]

iGSM_neu = (1/RGSM)*uGSM-(kGSM/RGSM)*wGSM;
fM = subs(fM, {iGSM}, iGSM_neu);

%% 2 Ruhelagen bestimmen, System linearisieren

%stationäre Eingangsgroeßen
uGSMR = 5.6;
MextR = 0;

% Ruhelagen bestimmen (Differentialgleichungen Null setzen)
eqns = [fM(1) == 0, fM(2) == 0, fM(3) == 0];
vars = [phiGSMP, wGSM, wP];
[phiGSMP_R, wGSM_R, wP_R]  = solve(eqns, vars)

rl_sym = subs([phiGSMP_R, wGSM_R, wP_R], {uGSM,Mext}, [uGSMR, MextR])
rl = double(subs(rl_sym, p_syms, p))

%Linearisieren um die Ruhelage
A_sym = jacobian(fM, xM);
% x_R und u_R einsetzen
A_sym = subs(A_sym,{phiGSMP, wGSM, wP},[phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
A_sym = subs(A_sym,{uGSM, Mext},[uGSMR, MextR])
A = double(subs(A_sym, p_syms, p))

%B-Eingangsvektor um 
Bu_sym = diff(fM, uGSM);
Bu_sym = subs(Bu_sym,{phiGSMP, wGSM, wP},[phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
Bu_sym = subs(Bu_sym,{uGSM, Mext},[uGSMR, MextR]);

%B-Eingangsvektor dm
Bd_sym = diff(fM, Mext);
Bd_sym = subs(Bd_sym,{phiGSMP, wGSM, wP},[phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
Bd_sym = subs(Bd_sym,{uGSM, Mext},[uGSMR, MextR]);

B = [double(subs(Bu_sym, p_syms, p)) double(subs(Bd_sym, p_syms, p))]

%C
C = [0,0, 0, 1]

dxM = A_sym*xM + Bu_sym*uGSM + Bd_sym*Mext;
y = C * wP;

%Eigenwerte
%Aeigen = double(eig(A))

%% 3. Ist die Ruhelage des reduzierten Systems gleich jener des vollstaendigen Systems?

% JA

% weil iGSM = (1/RGSM)*uGSM-(kGSM/RGSM)*wGSM genau der Ruhelage der DGL
% iGSM' = (1/LGSM)*(uGSM-RGSM*iGSM-kGSM*wGSM) entspricht

% iGSM' = (1/LGSM)*(uGSM-RGSM*iGSM-kGSM*wGSM) = 0
% uGSM - RGSM*iGSM - kGSM*wGSM) = 0
% uGSM - kGSM*wGSM = RGSM*iGSM
% (1/RGSM)*uGSM - (kGSM/RGSM)*wGSM = iGSM