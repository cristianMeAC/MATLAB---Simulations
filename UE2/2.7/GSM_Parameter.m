%%
%Andreas Stadelmann 01525173
%Cristian Avram     01304470

clear all   % Es werden alle Variablen gelöscht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zurückgesetzt

%% Systemparameter

LGSM = 1.4e-3;
RGSM = 0.46;
kGSM = 0.1;
JGSM = 12.4e-3;
dcGSM = 0.152;
dvGSM = 1.8e-3;
JP = 32.5e-3;
dcP = 0.169;
dvP = 2.7e-3;
dqP = 1e-4;
cGSMP = 0.6822;
dGSMP = 1e-5;

%% Anfangszustaende
% Verwenden Sie als Anfangszustand die Ruhelage für uGSM = 5.6V und Mext = 0Nm.
% Ruhelage aus 2.2.2

syms uGSM Mext iGSM phiGSMP wGSM wP

uGSMR = 5.6;
MextR = 0;

%Gleichungen 
MP = dcP + dvP*wP+dqP*wP^2+Mext;
MrGSM = dcGSM + dvGSM*wGSM;
%phiGSMP = phiGSM - phiP;
Mkopp = (wGSM-wP)*dGSMP + phiGSMP*cGSMP;
MGSM = kGSM*iGSM;

xM = [iGSM; phiGSMP; wGSM; wP];
fM = [(1/LGSM)*(uGSM-RGSM*iGSM-kGSM*wGSM); wGSM-wP; (MGSM-MrGSM-Mkopp)/(JGSM); (Mkopp-MP)/JP]

%Ruhelagen bestimmen
eqns = [fM(1) == 0, fM(2) == 0, fM(3) == 0, fM(4) == 0];
vars = [iGSM, phiGSMP, wGSM, wP];
[iGSM_R, phiGSMP_R, wGSM_R, wP_R]  = solve(eqns, vars)

ruhelage = double(subs([iGSM_R, phiGSMP_R, wGSM_R, wP_R], {uGSM,Mext}, [uGSMR, MextR]))

iGSM_0 = ruhelage(2,1);
phiGSMP_0 = ruhelage(2,2);
wGSM_0 = ruhelage(2,3);
wP_0 = ruhelage(2,4);

%% linearisiertes System

syms wP phiGSMP wGSM iGSM um dm

%Linearisieren
A = jacobian(fM, xM)
% x_R und u_R einsetzen
A = subs(A,{iGSM, phiGSMP, wGSM, wP},[iGSM_R(2), phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
A = subs(A,{uGSM, Mext},[uGSMR, MextR]);
A = double(A)

%B-Eingangsvektor um 
Bu = diff(fM, uGSM);
Bu = subs(Bu,{iGSM, phiGSMP, wGSM, wP},[iGSM_R(2), phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
Bu = subs(Bu,{uGSM, Mext},[uGSMR, MextR]);

%B-Eingangsvektor dm
Bd = diff(fM, Mext);
Bd = subs(Bd,{iGSM, phiGSMP, wGSM, wP},[iGSM_R(2), phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
Bd = subs(Bd,{uGSM, Mext},[uGSMR, MextR]);

B = [double(Bu) double(Bd)]

%C
C = [0,0, 0, 1]

%D
D = [0,0]

%%  linearisiertes reduziertes System

xM = [phiGSMP; wGSM; wP];
fM = [wGSM-wP; (MGSM-MrGSM-Mkopp)/(JGSM); (Mkopp-MP)/JP]

iGSM_neu = (1/RGSM)*uGSM-(kGSM/RGSM)*wGSM;
fM = subs(fM, {iGSM}, iGSM_neu);

%Ruhelagen bestimmen
eqns = [fM(1) == 0, fM(2) == 0, fM(3) == 0];
vars = [phiGSMP, wGSM, wP];
[phiGSMP_R, wGSM_R, wP_R]  = solve(eqns, vars)

ruhelage = double(subs([phiGSMP_R, wGSM_R, wP_R], {uGSM,Mext}, [uGSMR, MextR]))

%Linearisieren
A_red = jacobian(fM, xM)
% x_R und u_R einsetzen
A_red = subs(A_red,{phiGSMP, wGSM, wP},[phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
A_red = subs(A_red,{uGSM, Mext},[uGSMR, MextR]);
A_red = double(A_red)

%B-Eingangsvektor um 
Bu_red = diff(fM, uGSM);
Bu_red = subs(Bu_red,{phiGSMP, wGSM, wP},[phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
Bu_red = subs(Bu_red,{uGSM, Mext},[uGSMR, MextR]);

%B-Eingangsvektor dm
Bd_red = diff(fM, Mext);
Bd_red = subs(Bd_red,{phiGSMP, wGSM, wP},[phiGSMP_R(2), wGSM_R(2), wP_R(2)]);
Bd_red = subs(Bd_red,{uGSM, Mext},[uGSMR, MextR]);

B_red = [double(Bu_red) double(Bd_red)]

%C
C_red = [0, 0, 1]

%D
D_red = [0,0]


%% Regler

Sys = ss(A_red, B_red, C_red, D_red);
Gs = tf(Sys);

G = Gs(1)       %ohne Störung
Gd = Gs(2)      %mit Störung

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Regler entwerfen, Skriptum Seite 174

%(A)
%Vorgaben 
e = 0;          %die bleibende Regelabweichung e als Maß für die stationäre Genauigkeit
tr = 1;         %Anstiegszeit tr als Maß für die Schnelligkeit (Dynamik)
u = 0;          %die Überschwingweite M oder das prozentuelle Überschwingen u = (M-1)*100 als Maß für den Dämpfungsgrad (Dynamik)
Ta = 10e-3;     %Abtastzeit Ta = 10ms

%der offene Kreis muss mindestens eine einfache Polstelle bei s = 0 haben,
%damit e = 0

%(B)
%Abtastzeit Ta = 10ms
%w-Übertragungsfunktion G#(q) berechnen, s --> z --> q
%Transformation vom kontinuierlichen s-Bereich in den diskreten Zeitbereich
Gz = c2d(G, Ta, 'zoh')          %ohne Störung
Gdz = c2d (Gd, Ta, 'zoh')       %mit Störung
%Transformation vom diskreten s-Bereich in den kontinuierlichen q-Berecih
Gq = d2c(Gz, 'tustin')          %ohne Störung
Gdq = d2c(Gdz, 'tustin')        %mit Störung

%TI und xi bestimmen
%1+2xi(qT)+(qT)^2 das konjugiert komplexe Polpaar der Streckenübertragungsfunktion als Nullstellen besitzt
%und qr,j die gewünschten Realisierungspole bezeichnen

%zpkdata(SYS) returns the zeros, poles
[z,p] = zpkdata(Gq);
%die ersten 2 Pole
p1 = p{1}(1);
p2 = p{1}(2);

%aus den Polen 2 Nennerpolynome generieeren , (q-p1)(q-p2)
pn1 = [1 -p1];
pn2 = [1 -p2];

%Produkt bilden 
pr = conv(pn1, pn2);

%in Form 1 + 2*xi*T +(qT)^2 bringen --> durch die letzte Zahl dividieren 
pr = pr / pr(length(pr));

T = sqrt(pr(1));    % --> 0.115035427726282
xi = pr(2)/(2*T);   % --> 0.083270173687915

%(C)
%Formeln, Skriptum Seite 172, daraus folgt:
w0 = 2/Ta;    
wc = 1.2/tr;        %wc*tr =~ 1.2, wc..Durchtrittsfr. trennt jene Frequenzen die verstärkt werden, von jenen die abeschwächt werden (vom offenen Regelkreis) 
phi = 70-u;         %phi[%]+u[Grad] = 70, phi..Phasenreserve Maß für den Abstand zur Stabilitätsgrenze, Verminderung der Phasenreserve --> Zunahme der Schwingneigung bzw. des Überschwingens

%(D)
%R#q = (VI(1+q*TI)/q) * ((1+2*xi(q*T)+(q*T)^2)/(q-qr1)(q-qr2))

%Integrator hinzufuegen:
%Bodediagramm aller bekannten Terme des offenen Kreise
q = tf('s');
Rq1 = (1+2*xi*(q*T)+(q*T)^2)/q;
Lq1 = Rq1*Gq

%PI-Regler
RqPI1 = 1/q;
LqPI1 = RqPI1*Gq;
 
%aus (E)
%Um ein kriechendes Einlaufen der Sprungantwort in den stationären Endwert zu
%vermeiden, soll in (D) der Regler r#(q) so entworfen werden, dass ca. 1 Dekade
%um die Durchtrittsfrequenz wc die Betragskennlinie von L#(q) mit mindestens 20 dB/Dekade abfällt.
%dadurch folgt --> um Durchtrittsfrequenz herum +/-10grad
%atan(-wc/qr)=10grad, 
qr1 = -wc/tan(10*pi/180);
qr2 = qr1;
qrPr= (q-qr1)*(q-qr2);

Rq2 =(1+2*xi*(q*T)+(q*T)^2)/(q*qrPr);
Lq2 = Rq2*Gq

%An der Durchtrittsfrequenz wc errechnet sich das Argument von Lq2(Iwc) zu arg(Lq2(Iwc))
%arg(Lq2) = arctan(Im Zähler/Re Zähler) - arctan(Im Nenner/Re Nenner)
argLq2 = rad2deg(angle(evalfr(Lq2, i*wc)));
zn = argLq2 + 180;

%Damit muss wegen phi = 70 Grad mithilfe des Linearterms im Zähler
%des PI-Reglers (1+s*TI) die Phase um (70 - zn)  angehoben werden.
%daraus folgt --> arg(1+Iwc*TI) = arctan(wc*TI/1)-arctan(0)= zn*pi/180
an = phi - zn;
TI = (tan(an*(pi/180)))/(wc)

Rq3=(1+q*TI)*(1+2*xi*(q*T)+(q*T)^2)/(q*qrPr);
Lq3 = Rq3*Gq;

%PI-Regler
argLq2PI = rad2deg(angle(evalfr(LqPI1, i*wc)));
znPI = argLq2PI + 180;
anPI = phi - znPI;
TiPI= (tan(anPI*(pi/180)))/(wc)

RqPI2 = (1+q*TiPI)/q;
LqPI2 = RqPI2*Gq;

%In einem Schritt wird V so berechnet, dass der Amplitudengang
%bei wc die 0-dB-Linie schneidet,
%also die Bedingung VI|R(s)*G(s)|= 1 erfüllt.
VI = 1/abs(evalfr(Lq3, j*wc))

Rq4=VI*(1+q*TI)*(1+2*xi*(q*T)+(q*T)^2)/(q*qrPr);
Rq = Rq4;

%Regler als z-Übertragungsfunktion
Rz = c2d(Rq,Ta,'tustin');
[Rz_num, Rz_den] = tfdata(Rz);
Rz_num = Rz_num{1};
Rz_den = Rz_den{1};

%Regler als Differenzengleichung
RzSS=ss(Rz);
Ad=RzSS.a;
Bd=RzSS.b;
Cd=RzSS.c;
Dd=RzSS.d;

%PI-Regler
VIPI = 1/abs(evalfr(LqPI2, j*wc))

RqPI3 = VIPI*(1+q*TiPI)/q;
RqPI = RqPI3;

%Regler as z-Übertragungsfunktion
RzPI = c2d(RqPI,Ta,'tustin');
[RzPI_num, RzPI_den] = tfdata(RzPI);
RzPI_num = RzPI_num{1};
RzPI_den = RzPI_den{1};

%Regler als Differenzengleichung
RzPISS=ss(RzPI);
AdPI=RzPISS.a;
BdPI=RzPISS.b;
CdPI=RzPISS.c;
DdPI=RzPISS.d;