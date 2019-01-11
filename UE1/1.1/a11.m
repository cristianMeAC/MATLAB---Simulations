%%
%Andreas Stadelmann 01525173
%Cristian Avram     01304470

clear all   % Es werden alle Variablen gel�scht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zur�ckgesetzt
%%
%1.1.1. Gegeben ist das Gleichungssystem:
%    x1 + 2x2 + 4x3 = 5  (1.1a)
%   2x1 + 2x2 + x3  = 4  (1.1b)
%   3x1 + 2x2       = 1  (1.1c)
 
%Schreiben Sie dieses Gleichungssystem in Matrixdarstellung an und bestimmen
%Sie den L�sungsvektor x = [x1 x2 x3]T
% A * x = b 

A = [ 1 2 4; 2 2 1; 3 2 0];
b = [5; 4; 1];

%F�hren Sie die Rechnung einmal mit dem Befehl inv() und einmal mit dem
%Befehl mldivide() durch. Untersuchen Sie die Geschwindigkeits- und
%Genauigkeitsunterschiede der beiden Befehle. 

% Rechnung mit inv()
tic
x_inv = inv(A) * b
toc
norm_inv = norm(A * x_inv - b)

% Rechnung mit mldivide()
tic 
x_mldivide = A \ b
toc
norm_mldivide = norm(A * x_mldivide - b)

%Antwort:
%Die beiden Befehle sind auf meinem System meistens fast gleich schnell.
%Der numerische Fehler ist bei mldivide gleich 0. F�r das L�sen solcher 
%Systeme w�rde ich demnach den Befehl mldivide verwenden.

%%
%1.1.2. Gegeben ist die Fourier-Reihenentwicklung n-ter Ordnung mit
%Amplitude A = 10. Stellen Sie diese Rechteckfunktion f�r n = 1,2,...,100
%im Intervall x [0,10] mithilfe einer for-Schleife dar. Nutzen Sie zur
%schrittweisen grafischen Darstellung der Funktion in der for-Schleife den
%pause-Befehl. 

Amp = 10;
x_vec = 0:pi/1000:10;
rect_mat = [];

for n = 1:1:100
    
       % speichert den Summanden f�r k = n in rect_summand{n}
       rect_summand{n} = @(x) ( 4 / (pi*(2*n-1)) * sin((2*n-1)*x) );
       
       % erzeugt einen Vektor aus dem Funktion-Handle an{n} 
       rect_mat(n, :) = rect_summand{n}(x_vec);
       
       % summiert alle bisherigen Vektoren
       rect_vec = Amp * sum(rect_mat, 1);

       % plot Rechtecksfunktion
       figure(1)
       plot(x_vec, rect_vec);

       xlabel('x', 'FontName', 'Arial', 'FontSize', 14) 
       ylabel('Amplitude', 'FontName', 'Arial', 'FontSize', 14) 
      
       title(sprintf('Fourier Reihenentwicklung mit n=%d', n));

       % warte fuer 0.1s
       pause(0.1)
       
end
%%
% Aufgabe 1.3

t_konst = 0:0.1:2*pi;

x_vec = 0:pi/10:10;

% Koordinatenmatrix erzeugen
[xx,yy] = meshgrid(x_vec,x_vec);

f = @(x, y, t) sin(x*pi/10) .* sin(y*pi/10) * abs(sin(t));

for i = 1:1:length(t_konst)
    zz = f(xx, yy, t_konst(i));
    
    surf(xx, yy, zz);
    zlim manual
    zlim([0 1]);
    
    pause(0.1)
end



