%%
%Andreas Stadelmann 01525173
%Cristian Avram     01304470

clear all   % Es werden alle Variablen gel�scht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zur�ckgesetzt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%G(s) bestimmen
%Berechnungen aus a23
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

%die Ruhelage laut Angabe
um = 5.6;
dm = 0;

syms wGSM phiGSMP wP um dm
 
uGSM = um;
Mext = dm;

% Gleichungen
MP = dcP + dvP*wP+dqP*wP^2+Mext;
MrGSM = dcGSM + dvGSM*wGSM;
Mkopp = (wGSM-wP)*dGSMP + (phiGSMP)*cGSMP;
iGSM = (1/RGSM)*uGSM-(kGSM/RGSM)*wGSM;
MGSM = kGSM*iGSM;

fM_2 = [wGSM - wP; -(dcGSM + cGSMP*phiGSMP + dvGSM*wGSM + dGSMP*(wGSM - wP) - kGSM*(um/RGSM - (kGSM*wGSM)/RGSM))/JGSM; -(dcP + dm - cGSMP*phiGSMP + dvP*wP - dGSMP*(wGSM - wP) + dqP*wP^2)/JP];
xM = [phiGSMP; wGSM; wP];
phiGSMP =  0.506027297604110;
wGSM = 30.594990920230799;
wP = 30.594990920230799;

A = jacobian(fM_2, xM);
A = double(subs(A));

Bu = diff(fM_2, um);
Bu = subs(Bu);
Bd = diff(fM_2, dm);
Bd = subs(Bd);
B = [double(Bu) double(Bd)];

C = [0,0,1];

D= [0,0];

Sys = ss(A, B, C, D);
Gs = tf(Sys);

G = Gs(1)       %ohne St�rung
Gd = Gs(2)      %mit St�rung

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Regler entwerfen, Skriptum Seite 174

%(A)
%Vorgaben 
e = 0;          %die bleibende Regelabweichung e als Maß f�r die station�re Genauigkeit
tr = 1;         %Anstiegszeit tr als Ma� f�r die Schnelligkeit (Dynamik)
u = 0;          %die �berschwingweite M oder das prozentuelle �berschwingen u = (M-1)*100 als Ma� f�r den D�mpfungsgrad (Dynamik)
Ta = 10e-3;     %Abtastzeit Ta = 10ms

%der offene Kreis muss mindestens eine einfache Polstelle bei s = 0 haben,
%damit e = 0

%(B)
%Abtastzeit Ta = 10ms(laut Angabe)
%w-�bertragungsfunktion G#(q) berechnen, s --> z --> q
%Transformation vom kontinuierlichen s-Bereich in den diskreten Zeitbereich
Gz = c2d(G, Ta, 'zoh')          %ohne St�rung
Gdz = c2d (Gd, Ta, 'zoh')       %mit St�rung
%Transformation vom diskreten s-Bereich in den kontinuierlichen q-Berecih
Gq = d2c(Gz, 'tustin')          %ohne St�rung (von 2.3 die Strecke)
Gdq = d2c(Gdz, 'tustin')        %mit St�rung

%TI und xi bestimmen
%1+2xi(qT)+(qT)^2 das konjugiert komplexe Polpaar der Strecken�bertragungsfunktion als Nullstellen besitzt
%und qr,j die gew�nschten Realisierungspole bezeichnen

%zpkdata(SYS) returns the zeros, poles
[z,p] = zpkdata(Gq);
%die ersten 2 Pole
p1 = p{1}(1);  % -0.7239 + 8.6628i
p2 = p{1}(2);  % -0.7239 - 8.6628i 

%aus den Polen 2 Nennerpolynome generieeren , (q-p1)(q-p2)
pn1 = [1 -p1];
pn2 = [1 -p2];

%Produkt bilden(wir koennen nicht 2 Polynome multiplizieren, deshalb machen wir Faltung)
pr = conv(pn1, pn2);

%in Form 1 + 2*xi*T +(qT)^2 bringen --> durch die letzte Zahl dividieren 
pr = pr / pr(length(pr));

T = sqrt(pr(1));    % --> 0.115035427726282
xi = pr(2)/(2*T);   % --> 0.083270173687915

%(C)
%Formeln, Skriptum Seite 172, daraus folgt:
w0 = 2/Ta;    
wc = 1.2/tr;        %wc*tr =~ 1.2, wc..Durchtrittsfr. trennt jene Frequenzen die verst�rkt werden, von jenen die abeschw�cht werden (vom offenen Regelkreis) 
phi = 70-u;         %phi[%]+u[Grad] = 70, phi..Phasenreserve Ma� f�r den Abstand zur Stabilit�tsgrenze, Verminderung der Phasenreserve --> Zunahme der Schwingneigung bzw. des �berschwingens

%(D)
%R#q = (VI(1+q*TI)/q) * ((1+2*xi(q*T)+(q*T)^2)/(q-qr1)(q-qr2))

%Integrator hinzufuegen:
%Bodediagramm aller bekannten Terme des offenen Kreise
q = tf('s'); % damit MATLAB weiß dass wir q statt s in Uebertragungsfkt verwenden wollen
Rq1 = (1+2*xi*(q*T)+(q*T)^2)/q;
Lq1 = Rq1*Gq

%PI-Regler
RqPI1 = 1/q;
LqPI1 = RqPI1*Gq;
 
%aus (E)
%Um ein kriechendes Einlaufen der Sprungantwort in den station�ren Endwert zu
%vermeiden, soll in (D) der Regler r#(q) so entworfen werden, dass ca. 1 Dekade
%um die Durchtrittsfrequenz wc die Betragskennlinie von L#(q) mit mindestens 20 dB/Dekade abf�llt.
%dadurch folgt --> um Durchtrittsfrequenz herum +/-10grad
%atan(-wc/qr)=10grad, 
qr1 = -wc/tan(10*pi/180);
qr12 = -1000;
qr2 = qr1;
qrPr= (q-qr1)*(q-qr2);
qrPr2 = (q-qr12)*(q-qr12);

Rq2 =(1+2*xi*(q*T)+(q*T)^2)/(q*qrPr);
Lq2 = Rq2*Gq

figure
bode((1+2*xi*(q*T)+(q*T)^2)/(qrPr));
hold on
bode((1+2*xi*(q*T)+(q*T)^2)/(qrPr2));

%An der Durchtrittsfrequenz wc errechnet sich das Argument von Lq2(Iwc) zu arg(Lq2(Iwc))
%arg(Lq2) = arctan(Im Z�hler/Re Z�hler) - arctan(Im Nenner/Re Nenner)
argLq2 = rad2deg(angle(evalfr(Lq2, i*wc)));
zn = argLq2 + 180;

%Damit muss wegen phi = 70 Grad mithilfe des Linearterms im Z�hler
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
%also die Bedingung VI|R(s)*G(s)|= 1 erf�llt.
VI = 1/abs(evalfr(Lq3, j*wc))

Rq4=VI*(1+q*TI)*(1+2*xi*(q*T)+(q*T)^2)/(q*qrPr);
Lq = Rq4*Gq;

%PI-Regler
VIPI = 1/abs(evalfr(LqPI2, j*wc))

RqPI3 = VIPI*(1+q*TiPI)/q;
LqPI = RqPI3*Gq;

%Sprungantwort des geschlossnen Kreises
%L = zl/nl, geschlossener Kreis Try = L / (1+L) --> Try = zl/(zl+nl)
[num,den]=tfdata(Lq);
den = cell2mat(den);
num = cell2mat(num);
Tqry = tf (num, num + den)      %F�hrungs�bertragungssystem
Tqdy = Gdq/(1+Lq)               %St�r�bertragungssystem

%%%%%%%%%%%%
%PI-Regler
[numPI,denPI]=tfdata(LqPI);
denPI = cell2mat(denPI);
numPI = cell2mat(numPI);
TqryPI = tf (numPI, numPI + denPI)      %F�hrungs�bertragungssystem
TqdyPI = Gdq/(1+LqPI)               %St�r�bertragungssystem
%%
%Grafiken

figure;
bode(Lq1);
hold on
bode(Lq2);
bode(Lq3);
bode(Lq);
legend("Lq1", "Lq2", "Lq3", "Lq");

%Bodediagramm, PI-Regler vs. Kompensationsregler
figure;
bode(Lq);
hold on;
grid;
bode(LqPI);
title("PI-Regler vs. Kompensationsregler");
legend("L#(q)","L#(q) (PI)");

%F�hrungs�bertragungssystem vs. St�r�bertragungssystem
figure
step(Tqry)
hold on
grid on
step(TqryPI)


step(Tqdy)
step(TqdyPI)
title("F�hrungs�bertragungssystem, St�r�bertragungssystem"')
legend('Kompensationsregler - Tqry', 'PI-Regler - Tqry(PI)', 'Kompensationsregler - Tqdy', 'PI-Regler - Tqdy(PI)')
legend('Location','southeast')
