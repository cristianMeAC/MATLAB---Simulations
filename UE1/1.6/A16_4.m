%%
%Aufgabe 1.6.4 Datenmenge Sr f�r ur = 1 berechnen
clear all
close all
clc

Ta = 0.000001;
T1 = 0.001;
c = 0.5;
a = 0.4;

phi = [exp(-Ta/T1), c*a*(1-exp(-Ta/T1)); Ta, 1-a*Ta];
E = [1,0;0,1];
gamma = [0;Ta];

inv(E-phi)*gamma