%%
% Aufagbe 1.6.3 
% Prozessor Parameter
clear all
close all 
clc

T1 = 0.001;
Ta = 0.000001;
c = 0.5;
a = 0.4;

% Zustandsraumdarstellung Prozessor 
A = [exp(-Ta/T1)];
B = [(1-exp(-Ta/T1))*c];
C = [1];
D = 0;

%Aufgabe 1.6.5 

Ts = 15e-3;
Te = 10e-3;

Phi = [exp(-Ta/T1),c*a*(1-exp(-Ta/T1));Ta,1-a*Ta];
Gamma = [0;Ta];

k = Ts/Ta;
Sr = [1;5];
Phi_k = Phi^k

t=0:Ta:Ts;
u=heaviside(t-Te)*1000;
plot(t,u)

sum = 0;
for j = 0:1:k-1
    sum = sum + Phi^(k-j-1) * Gamma * u(j+1);
end

Sk = Phi_k * Sr + sum