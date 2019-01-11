% Parameterfiles fuer das Simulink Modell des elektrischen Systems.
close all;
clear all;
clc;

% Widerstaende R1 und R2
R_1 = 625;
R_2 = 3500;

% Widerstandsverhaeltnis R4
K = 3;

% Kapazitaeten C1 und C2
C1_ref  = 1e-6;
UC1_ref = -10;
kC1     = 800e-9;
C2      = 1e-6;

% Festlegung der Ruhelage
UeR = 5;
UsR = 5;

%% Eingangsspannung und Stoerspannung
% Sprungfoermige Eingangsspannung
te = .05;      % Einschaltzeitpunkt
Ue = 7;      % Endwert des Sprunges

% Sinusfoermige Eingangsspannung
Usinus = 1;    % Amplitude des Sinus    
we     = 1000; % Winkelfrequenz des Sinus   

% Sprungfoermige Stoerspannung
ts = .1;     % Einschaltzeitpunkt   
Us = 4;     % Endwert des Sprunges
% Rampenförmige Störspannung
tsr = 0.15; %Einschaltzeitpunkt
dr = 10;    % Steigung


%% Modell des linearisierten Systems

% Aufgabe 1.4.4
% -------------
% Errechnen Sie die Ruhelage des Systems!
% UC1R = ;
% UC2R = ;
% UaR  = ;

% Errechnen Sie für die gegebene Ruhelage die Ableitungen der in C1 
% gespeicherter Ladung Q1 nach der anliegenden Spannung!
% dQ1  = ; % Erste Ableitung von Q1 nach UC1
% ddQ1 = ; % Zweite Ableitung von Q1 nach UC1

% Ergaenzen Sie die Systemmatrix A (Asys), die Eingangsvektoren bu (busys)
% und bd (bdsys) sowie den Ausgangsvektor c (csys) und den Durchgriff d (dsys)
% für das linearisierte System!

% A11 = ;
% A12 = ;
% A21 = ;
% A22 = ;

% A  = [A11,A12;A21,A22];
% bu = [];
% bd = [];
% c  = [];
% du = [];
% dd = [];

