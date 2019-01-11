%%
%Andreas Stadelmann 01525173
%Cristian Avram     01304470

clear all   % Es werden alle Variablen gel�scht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zur�ckgesetzt

%zum Testen von verschiedenen Anstiegszeiten die Variable tr in GSM_Parameter.m �ndern
run("GSM_Parameter.m");

% �bertragungsfunktion vom Sensorrasuchen n zum Eingang u:
Tnu = (-Rq) / (1+Gq*Rq);

% �bertragungsfunktion vom Sensorrasuchen n zum Eingang u:
Tny = (-Rq*Gq) / (1+Gq*Rq);

figure
bode(Tny)
hold on
bode(Tnu)
legend("Tny(q)", "Tnu(q)")

% bei kleinen Anstiegszeiten (tr=0.1) ist deutlich zu erkennen, dass die Stellgr��e durch die
% Verst�rkung von Tnu im hochfrequenten stark beeinflusst wird
