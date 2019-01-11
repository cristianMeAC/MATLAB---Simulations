% Parameterfiles fuer das Simulink Modell des elektrischen Systems 
% (Aufgabe 1.4 und 1.5)
close all;
clear all;
clc;

% Parameter
R1      = 625;    % Widerstand in Ohm
R2      = 3500;   % Widerstand in Ohm
K       = 3;      % Spannungsverst‰rkung
C1_ref  = 1e-6;   % Referenz-Kapazit‰tswert auf Kennlinie von C1 in F
UC1_ref = -10;    % Referenz-Spannungswert auf Kennlinie von C1 in V
kC1     = 800e-9; % Steigung der Kennlinie von C1 in F/V
C2      = 1e-6;   % Kapazit‰t in F

% Festlegung der Ruhelage in Volt
UeR = 5;
UsR = 5;

%% Eingangsspannung und Stoerspannung (Aufgabe 1.5)
% Sprungfoermige Eingangsspannung
te = 0.05;      % Einschaltzeitpunkt in Sekunden
Ue = 7;      % Endwert des Sprunges in Volt

% Sinusfoermige Eingangsspannung
Usinus = 1; % Amplitude des Sinus in Volt
we     = 1000;  % Winkelfrequenz des Sinus in rad/s

% Sprungfoermige Stoerspannung
ts = 0.1;     % Einschaltzeitpunkt in Sekunden
Us = 4;     % Endwert des Sprunges in Volt

% Rampenfˆrmige Stˆrspannung
tsr = 0.15; %Einschaltzeitpunkt
dr = 10;    % Steigung

%% Modell des linearisierten Systems

% Aufgabe 1.4.4
% -------------
% Errechnen Sie die Ruhelage des Systems!
UC1R = -(K-2)*(R2*UeR+(R1+R2)*UsR)/(R1+2*R2)-UsR;
UC2R = (R2*UeR+(R1+R2)*UsR)/(R1+2*R2);
UaR  = K * UC2R;

% Errechnen Sie f¸r die gegebene Ruhelage die Ableitungen der in C1 
% gespeicherter Ladung Q1 nach der anliegenden Spannung!
dQ1  = C1_ref-kC1*UC1_ref+kC1*UC1R; % Erste Ableitung von Q1 nach UC1
ddQ1 = kC1; % Zweite Ableitung von Q1 nach UC1

% Ergaenzen Sie die Systemmatrix A (Asys), die Eingangsvektoren bu (busys)
% und bd (bdsys) sowie den Ausgangsvektor c (csys) und den Durchgriff d (dsys)
% f¸r das linearisierte System!

dUC1 = 1/dQ1 * (UeR/R1 - UC1R*(R1+R2)/(R1*R2) - UC2R*(K*(R1+R2)-R1)/(R1*R2));

A11 = 1/dQ1 * (ddQ1*dUC1 - (R1+R2)/(R1*R2));
A12 = -(K*(R1+R2)-R1)/(dQ1*R1*R2);
A21 = 1/(C2*R2);
A22 = (K-2)/(C2*R2);

A  = [A11,A12;A21,A22];
bu = [1/(dQ1*R1);0];
bd = [0;1/(C2*R2)];
B = [bu,bd];
c  = [0,K];
du = 0;
dd = 0;
%% Uebertragungsfunktion G (Eingang u -> Ausgang y) und Uebertragungsfunktion Gd (Stoerung d -> Ausgang y)

% Aufgabe 1.4.4
% -------------
% Bestimmen Sie zun‰chst eine MISO- oder zwei SISO-Zustandsraum- 
% darstellungen mittels ss(). Anschlieﬂend koennen Sie in beiden Faellen
% die gesuchten ‹bertragungsfunktionen (G und Gd) mittels tf() bestimmen.

sys_lin = ss(A,B,c,0);
G2 = tf(sys_lin);

G  = G2(1)
Gd = G2(2);

%% Verstaerkungsfaktor V, Daempfungsgrad xi und Zeitkonstante T von G

% Aufgabe 1.4.5
% -------------
% Bestimmen Sie den Verstaerkungsfaktor V, den Daempfungsgrad xi und die
% Zeitkonstante T der ‹bertragungsfunktion G. 

V  = (1.371e06/9.959e05);
T  = sqrt(1/9.959e05);
xi = (1600)/(9.959e05*2*T);

%% Bodediagramme der Uebertragungsfunktionen G und Gd

% Aufgabe 1.4.6
% -------------
% Zeichnen und interpretieren Sie die Bodediagramme der beiden
% Uebertragungsfunktionen G und Gd. Verwenden Sie dazu den Befehl bode().

figure 
bode(G)
hold on
bode(Gd)
hold on 
legend("G(s)","Gd(s)");