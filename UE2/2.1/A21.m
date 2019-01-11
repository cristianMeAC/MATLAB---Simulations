%%
%Andreas Stadelmann 01525173
%Cristian Avram     01304470

clear all   % Es werden alle Variablen gel�scht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zur�ckgesetzt

run("Esys_Parameter.m");

%% 1. Pi-Regler

% A 
%Vorgaben

V  = 1.377;
T  = 0.001002;
xi = 0.802;
G = tf ([V],[T^2 2*xi*T 1]);     %�bertragungsfunktion

% B 
%Vorgaben an den Frequenzgang des offenen Kreises

tr = 0.003;
ue = 5;
e = 0;

%Daraus folgt:

%wc*tr = 1.5 rad/s und %phi[Grad]+ue[%]=70 (Skriptum Seite 121) -->  
wc = 1.5/tr;     % 500rad/s
phi = 70 - ue;    % 65Grad 

%Damit die Regelabweichung e fuer r(t) = sigma(t) Null ist,
%muss der offene Regelkreis einen Pol bei s = 0 haben.

% C
%Als Regler in der Form R(s) = (VI*(1+s*TI))/s
%Im ersten Schritt: Bodediagramm aller bekannten Terme des offenen Kreises
%Ls = R(s)*G(s)

L1 = tf ([V], [T^2 2*xi*T 1 0]);

%An der Durchtrittsfrequenz wc = 500rad/s errechnet sich das Argument
%von L1(Iwc) zu arg(L1(Iwc))
%arg(L1) = arctan(Im Zaehler/Re Zaehler) - arctan(Im Nenner/Re Nenner)

argL1 = rad2deg(angle(evalfr(L1, i*wc)));
zn = argL1 + 180;

%Damit muss wegen phi = 65Grad mithilfe des Linearterms im Zaehler
%des PI-Reglers (1+s*TI) die Phase um (65Grad - zn))  angehoben werden.
%Daraus folgt 
%   arg(1+Iwc*TI) = arctan(wc*TI/1)-arctan(0)= zn*pi/180

an = phi - zn;
TI = (tan(an*(pi/180)))/(wc);

%In einem Schritt wird V so berechnet, dass der Amplitudengang
%bei wc die 0-dB-Linie schneidet,
%also die Bedingung VI|R(s)*G(s)|= 1 erfuellt.

L2 = tf ([V*TI V], [T^2 2*xi*T 1 0]);
VI = 1/abs(evalfr(L2, j*wc));

R_pi = tf ([VI*TI VI], [1 0]);

%-->daraus folgt L(s) = R(s)*G(s)

L_pi = series(R_pi, G);
L = tf ([V*VI*TI V*VI], [T^2 2*xi*T 1 0]);

%Sprungantwort des geschlossnen Kreises
%L = zl/nl, geschlossener Kreis Try = L / (1+L)
%--> Try = zl/(zl+nl)

% [num,den]=tfdata(L);
% den = cell2mat(den);
% num = cell2mat(num);
% Try = tf (num, num + den);
Try = feedback(L_pi, 1);

%Grafiken
opts=bodeoptions('cstprefs');
opts.XLim={[10 10^4]};
opts.grid = 'on';

figure
bode(L_pi,opts)
hold on
bode(L1,opts)
legend("L(s)","L1(s)");

figure 
step(Try)

% D
%Aus dem Bodediagramm von L(s) erkennt man,
%dass die Bedingungen von (B) erfuellt sind.
%L(s) faellt um wc mit mindestens 20dB/Dekade.

% E 
%Die zugehoerige Sprungantwort des geschlossenen Kreises zeigt,
%dass die Entwurfsanforderungen recht gut erfuellt werden.

% Regelabweichung --> 0 , ueberschwingungen beginnen bei tr
% u = (M-1)*100 --> Maximum = (u+100)/100 = 1.05
max = (ue+100)/100


% F 
%In diesem Beispiel wurden an die Stellgroe�e keine Anforderungen gestellt.


%% 2. Ist das Frequenzkennlinienverfahren ein exaktes Entwurfsverfahren? 

% Nein, wenn die Fuehrungsuebertragungfunktion ein PT2 Glied ist, dann Ja.
% Bei uns ist jedoch die Strecke ein PT2 Glied, also Nein. 

%% 3. PID-Regler

% Uebertragungsfunktion eines PID-Reglers  (S. 86)
% R(s) = Vp*(1+Ti*s)*(1+Td*s) / s*(1+Tr*s)

% weitere Anforderungen
e_pid = 1*10^-3;

%in der Angabe
Ti_pid = 1.1 * 10^-3;

% D-T1-Glied verhaelt sich annaehernd wie D-Glied fuer w << 1/Tr (Seite 84)
  
% mind. doppelte Nullstelle bei s=0 damit die Regelabweichung bei r(t) = t Null wird (Seite 123)
% bei einfacher Nullstelle ist e_inf = 1/V bei r(t) = t (PID hat einfache
% Nullstelle)
% daraus folgt Vl >= 1/e_inf damit die Regelabweichung eingehalten wird.
% (Seite 127)
Vl_pid = 1/e_pid;
% Damit errechnet sich unmittelbar der Verstaerkungsfaktor des Reglers zu
Vp_pid = Vl_pid/V

% Betrag und Phase von L1 = R1 * G an der Stelle w = wc betrachten
% offener Kreises anhand des bekannten PI-Reglers
R1_pid = tf([Ti_pid*Vp_pid Vp_pid], [1 0])
L1_pid = series(R1_pid, G);

%figure
%bode(L1_pid,opts)
%hold on
%bode(L, opts)
%legend("L1-pid", "L");

L1_pid_fr = evalfr(L1_pid, 1i*wc);
absL1_pid_wc = abs(L1_pid_fr)
argL1_pid_wc = angle(L1_pid_fr);
radtodeg(argL1_pid_wc)

% Phase beim wc ist negativ => Lag Glied(Wiki)
% Betrag und Phase muessen gesenkt werden -> Lag-Glied
% R_lag = (1 + s*T_lag) / (1 + s*eta*T_lag), eta > 1 (S. 128)

da = 1/absL1_pid_wc;
dphi = (-180 + phi) - radtodeg(argL1_pid_wc);

% Formeln zur Berechnung von T_lag und eta stehen im Skriptum auf S. 129
% Td = T_lag, Tr = eta*T_lag
T_lag = (da * sqrt(1+tand(dphi)^2) -1) / (wc * tand(dphi))
eta = (wc * T_lag - tan(dphi)) / (wc*T_lag *( 1 + wc*T_lag*tand(dphi)))
% zu ungenau -> fsolve
% fsolve solves a problem specified by F(x) = 0 for x, where F(x) is a function 
% that returns a vector value.
% Bedingungen stehen auf S. 129
% @solveLag - eigene Funktion
x_lag = fsolve(@solveLag, [T_lag,eta]);
T_lag = x_lag(1);
eta = x_lag(2);

Td_pid = T_lag;
Tr_pid = T_lag*eta;

% uebertragungsfunktionen erzeugen und darstellen 
R_lag = tf([Td_pid 1], [Tr_pid 1])
R_pid = series(R1_pid, R_lag);
L_pid = series(R_pid, G)

figure
bode(L_pi,opts)
hold on
bode(L_pid,opts)
hold on
bode(L1,opts)
legend("L(s)", "Lpid(s)", "L1(s)");

%{
% Td so bestimmen damit Phase um an_pid angehoben wird
L1_fr = evalfr(L1,1i*wc);
argL1_pid = angle(L1_fr);
an_pid = phi + (-180 - radtodeg(argL1_pid));
argTi = angle(1+i*wc*Ti_pid);
Td_pid = (tan((an_pid)*(pi/180)))/(wc);

% arg((1+Ti*iw)*(1+Td*iw)) = arg((1-Ti*Td*w^2) + (Ti*w+Td*w)*i)
% = arctan( (Ti*w+Td*w) / (1- Ti*Td*w^2) ) = an * pi/180
tan_an = tan(an_pid*(pi/180));
%Td_pid = (tan_an - Ti_pid * wc) / (wc + Ti_pid * wc^2 * tan_an);

% V so berechnet, dass der Amplitudengang bei wc die 0-dB-Linie schneidet
zl_pid = V * (1+Ti_pid*i*wc) * (1+Td_pid*i*wc);
nl_pid = (1+2*xi*i*wc*T + (i*wc*T)^2) * i*wc *(1+Tr_pid*i*wc)

Vp_pid = abs(nl_pid/zl_pid);



R_pid = tf ([Ti_pid*Td_pid Ti_pid+Td_pid 1], [Tr_pid 1 0]);
L_pid = R_pid * G;
Gg_pid = feedback(L_pid, 1);

figure
bode(L,opts)
hold on
bode(L_pid,opts)
hold on
bode(L1,opts)
legend("L(s)", "Lpid(s)", "L1(s)");

figure
step(Gg_pid);
grid on

%}

%%

% Logic 0 is not stable
isstable(L_pid)

%(Nyquist-Kriterium in Frequenzkennliniendarstellung) (S.171)
%A.. Der Verstaerkungsfaktor V ist positiv
% --> V = 1.377, efuellt
V

%B.. grad(nL(s)) + roh > grad(zL(s))
% --> 4+2 > 2, erfuellt
L_pid
% 
% %C.. nL und zL Hurwitzpolynome, und roh aus 0,1,2
% % in Hurwitzpolynom ist ein reelles Polynom, dessen Nullstellen alle einen echt negativen Realteil haben.
% % --> roh = 2, erfuellt
% % -->   -8.0040 + 5.9613i
% %       -8.0040 - 5.9613i
% %       -0.5454 + 0.0000i, erfuellt
% % [NUM,DEN] = tfdata(SYS) returns the numerator(s) and denominator(s) of the transfer function SYS
% [num,den]=tfdata(L_pid);
% %nenner 
% den = cell2mat(den);
% syms s
% hurwitz=poly2sym(den(1:end-1),s)
% hurroots=roots(den(1:end-1))
% %zaehler
% num = cell2mat(num);
% hurwitz=poly2sym(num(1:end-1),s)
% hurroots=roots(num(1:end-1))
% 
% 
% %D.. die Betragskennlinie von L(Iw) weisst genau einen Schnittpunkt mit der 0-dBLinie (eine Durchtrittsfrequenz w) auf
% figure
% bodemag(L_pid)
% hold on
% grid on
% plot(wc,0,'*')
% 
% %E.. im Bereich |L(Iw)|dB = 0 gelte -540 Grad < arg(L(Iw)) < 180 Grad (d. h. die Ortskurve
% %des offenen Kreises L(s) kann vor ihrem Eintauchen in den Einheitskreis den
% %Nullpunkt hoechstens einmal vollstaendig umkreisen).
% % erfuellt
% figure
% nyquist(L_pid)
% hold on
% 
% % --> alle Kriterien erfuellt
% %Unter diesen Voraussetzungen ist der Regelkreis mit der Uebertragungsfunktion
% %des offenen Kreises L(s) genau dann BIBO-stabil, wenn der Abstand der
% %Phase an der Durchtrittsfrequenz arg(L(Iwc)) zu -pi, die sogenannte Phasenreserve phi,
% %phi = arg(L(Iwc)) + pi
% %positiv ist.
% 
% phi = angle(evalfr(L1_pid, 1i*wc)) + pi

%% 4. Fuehrungsuebertragungsfunktion, Stoeruebertragungsfunktionen und Sprungantworten

% Fuehrungsuebertragungsfunktion des geschlossenen Kreises: Try = L / (1+L)
Try_pi = feedback(L_pi, 1);
Try_pid = feedback(L_pid, 1);

% Stoeruebertragungsfunktionen Tdy = Gd / (1+L)
Tdy_pi = Gd / (1+L_pi);
Tdy_pid = Gd / (1+L_pid);

%Sprungantworten

figure
step(Try_pi);
hold on
step(Try_pid);
grid on
legend("Try-pi", "Try-pid");

figure
step(Tdy_pi);
hold on
step(Tdy_pid);
grid on
legend("Tdy-pi", "Tdy-pid");


