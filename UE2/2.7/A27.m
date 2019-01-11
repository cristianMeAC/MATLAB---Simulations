%%
%Andreas Stadelmann 01525173
%Cristian Avram     01304470

clear all   % Es werden alle Variablen gelöscht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zurückgesetzt

%zum Testen von verschiedenen Anstiegszeiten die Variable tr in GSM_Parameter.m ändern
run("GSM_Parameter.m");

% Übertragungsfunktion vom Sensorrasuchen n zum Eingang u:
Tnu = (-Rq) / (1+Gq*Rq);

% Übertragungsfunktion vom Sensorrasuchen n zum Eingang u:
Tny = (-Rq*Gq) / (1+Gq*Rq);

figure
bode(Tny)
hold on
bode(Tnu)
legend("Tny(q)", "Tnu(q)")

% bei kleinen Anstiegszeiten (tr=0.1) ist deutlich zu erkennen, dass die Stellgröße durch die
% Verstärkung von Tnu im hochfrequenten stark beeinflusst wird
