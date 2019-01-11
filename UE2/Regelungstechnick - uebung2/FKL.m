%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aufgabe 2.1
% FKL Reglerentwurf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anforderungen an den offenen und geschlossenen Kreis
tr=3e-3;
ue=5;


%% Aufgabe 2.1.2 PI Regler

%Gesamtregler und offener Kreis
% R_PI = 
% L_PI = 
% bode(L_PI,'r');

%% Aufgabe 2.1.3 PID Regler

% Bleibende Regelabweichung
e_inf = 1e-3;

% Zeitkonstante des Intagralterms
TI = 1.1e-3;

%Gesamtregler und offener Kreis
% R_PID = 
% L_PID = 
% bode(L_PID,'g');

%% Aufgabe 2.1.4
% Kontrolle des Verhaltens des geschlossenen Kreise

T = 0.05;   % Simulationsdauer
figure;
subplot(3,1,1);hold on;grid on;
title('Sprungantwort Führungsübertragungsfunktion')

subplot(3,1,2);hold on;grid on;
title('Sprungantwort Störübertragungsfunktion')

subplot(3,1,3);hold on;grid on;
title('Rampenantwort Störübertragungsfunktion')