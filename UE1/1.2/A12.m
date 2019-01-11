%%
%Andreas Stadelmann 01525173
%Cristian Avram 01304470

clear all   % Es werden alle Variablen gel�scht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zur�ckgesetzt

%1.2 Betrachten Sie die Differentialgleichung zweiter Ordnung zur
%Beschreibung eines hamonischen Oszillators. Schreiben Sie ein m-file, in
%welchem der Integrationsalgorithmus auf das Differentialgleichungssystem
%angewendet wird. Variieren Sie die Abtastzeit und untersuchen Sie deren
%Einfluss auf die Genauigkeit der numerisch berechneten L�sung. 

%Anfangsbedingung
x0 = [0;0];
%Dynamikmatrix A 
A = [0,1;-2,-1];
%Eingangsmatrix B
B = [0;1];

%explizites Euler-Verfahren mit Abtastzeit Ta = 0.05s 
Ta1 = 0.05;
range1 = 0:Ta1:20;
u = ones(length(range1),1);  %Einheitssprung u
x1 = x0(:,1);
for k=2:length(range1)      
    x1(:,k) = x1(:,k-1)+Ta1*(A*x1(:,k-1)+B*u(k-1));
end
euler_steps1 = (x1(1,:));     %x(1,:) gibt die erste Reihe 

%explizites Euler-Verfahren mit Abtastzeit Ta = 0.01s 
Ta2 = 0.01;
range2 = 0:Ta2:20;
u = ones(length(range2),1);  %Einheitssprung u
x2 = x0(:,1);
for k=2:length(range2)      
    x2(:,k) = x2(:,k-1)+Ta2*(A*x2(:,k-1)+B*u(k-1));
end
euler_steps2 = (x2(1,:));     %x(1,:) gibt die erste Reihe 

%explizites Euler-Verfahren mit Abtastzeit Ta = 0.03s 
Ta3 = 0.03;
range3 = 0:Ta3:20;
u = ones(length(range3),1);  %Einheitssprung u
x3 = x0(:,1);
for k=2:length(range3)      
    x3(:,k) = x3(:,k-1)+Ta3*(A*x3(:,k-1)+B*u(k-1));
end
euler_steps3 = (x3(1,:));     %x(1,:) gibt die erste Reihe 

%explizites Euler-Verfahren mit Abtastzeit Ta = 0.1s 
Ta4 = 0.1;
range4 = 0:Ta4:20;
u = ones(length(range4),1);  %Einheitssprung u
x4 = x0(:,1);
for k=2:length(range4)      
    x4(:,k) = x4(:,k-1)+Ta4*(A*x4(:,k-1)+B*u(k-1));
end
euler_steps4 = (x4(1,:));     %x(1,:) gibt die erste Reihe 

%explizites Euler-Verfahren mit Abtastzeit Ta = 0.001s 
Ta5 = 0.001;
range5 = 0:Ta5:20;
u = ones(length(range5),1);  %Einheitssprung u
x5 = x0(:,1);
for k=2:length(range5)      
    x5(:,k) = x5(:,k-1)+Ta5*(A*x5(:,k-1)+B*u(k-1));
end
euler_steps5 = (x5(1,:));     %x(1,:) gibt die erste Reihe 

%Control System Toolbox
C = [1,0];
%Durchschaltmatrix
D = 0;

%generiert Kontiunierliches System
sys = ss(A,B,C,D);
%System von kontinuierlich zu diskret, c2s(sys,Ts)
sysd = c2d(sys,Ta1);
%Sprungantwort, step(dynamic sys, TFinal)
stepd = step(sysd,range1(end)); 

%Grafik
figure 
hold on
plot(range1,euler_steps1,'k.')
hold on 
plot(range2,euler_steps2,'b.')
hold on 
plot(range3,euler_steps3,'g.')
hold on 
plot(range4,euler_steps4,'r.')
hold on 
plot(range5,euler_steps5,'c.')
hold on
stepplot(sys,range1(end))

leg_range1 = 'Abtastzeit Ta = 0.05s';
leg_range2 = 'Abtastzeit Ta = 0.01s';
leg_range3 = 'Abtastzeit Ta = 0.03s';
leg_range4 = 'Abtastzeit Ta = 0.1s';
leg_range5 = 'Abtastzeit Ta = 0.001s';
leg_sys = 'Control System Toolbox';

legend(leg_range1,leg_range2,leg_range3,leg_range4,leg_range5,leg_sys)
hold off