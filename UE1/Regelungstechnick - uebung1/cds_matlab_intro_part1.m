%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Einf�hrung in Matlab/Simulink
% 
%
% Institut f�r Automatisierungstechnik 
% Gruppe f�r komplexe dynamische Systeme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dieses Musterfile dient zum Erlernen der grundlegenden Matlabbefehle. Es
% handelt sich dabei um ein sogenanntes m-file, in welchem der Code zu
% einer Datei zusammengefasst und direkt ausgef�hrt werden kann. Neben
% diesem m-file kann der Code auch direkt im Command Window hinter der
% Eingabeaufforderung ">>" eingegeben und ausgef�hrt werden. 

% Das Prozentzeichen dient in Matlab zum Einbinden von Kommentaren.
% Weiterhin ist es m�glich durch zwei aufeinanderfolgende Prozentzeichen
% den Quelltext in Abschnitte zu unterteilen. Das Ausf�hren der Datei
% erfolgt entweder durch Eingabe des Dateinamens im Command Window oder
% durch Dr�cken der F5-Taste, einzelne Abschnitte k�nnen durch die
% Kombination STRG+Eingabe ausgef�hrt werden. Einzelne Zeilen im
% Programmcode k�nnen weiterhin markiert und durch Dr�cken der F9-Taste
% ausgef�hrt werden.
%
% Das Ergebnis der Auswertungen erscheint im Command Window. Es wird dabei
% zeilenweise ausgewertet. Wird die Programmzeile mit einem Semikolon
% abgeschlossen, so wird die Ausgabe unterdr�ckt.

%%

% Beim Start eines neuen Programmes sollten zun�chst einmal alle alten
% Variablen gel�scht und alle offenen Fenster gel�scht werden

clear all   % Es werden alle Variablen gel�scht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zur�ckgesetzt

%%

% Matlab Hilfe: Es gibt zu allen Matlabbefehlen eine Hilfe, die entweder
% direkt im Command Window oder aber in einem zus�tzlichen Hilfefenster
% angezeigt wird (die Hilfe enth�lt des Weiteren auch sehr gute
% Einf�hrungsbeispiele

help sqrt % Hilfe im Command Window
%doc sqrt  % Hilfe in der Dokumentation (meist ausf�hrlicher)

%%

% Zuweisungen in Matlab: Matlab verwendet standardm��ig Flie�kommazahlen
% (double: 64 Bit).

variable_1 = 2;    % Durch das Semikolon wird die Ausgabe unterdr�ckt
variable_2 = 2.0;  % Kein Unterschied zu vorher
variable_3 = 4e5   % Durch Weglassen des Semikolons wird 
                   % das Ergebnis der Zeile ausgegeben
variable_3  = variable_1 % �berschreiben von Variablen                 
                   
% Einige Variablen sind vordefiniert

variable_4 = pi    % Kreiszahl pi
variable_5 = i     % Komplexe Einheit
variable_6 = j     % Komplexe Einheit
variable_7 = eps   % Flie�kommagenauigkeit
variable_8 = NaN   % Not a number, ung�ltiges Ergebnis 

% Zuweisen von String-Variablen

variable_9 = 'Hello World'  % Eingabe von Namen mit Hochkomma
variable_10 = strcat('Hello', ' World') % Zusammenf�gen von Strings

%%

% Mathematische Funktionen und Operatoren 

v1 = 2;
v3 = 3;

v1+v3   % Addition
v1*v3   % Multiplikation
v1-v3   % Subtraktion
v1/v3   % Division

mod(v1,v3)  % Modulodivision
rem(v1,v3)  % Rest der Division

sqrt(v1)        % Wurzel
exp(v1)         % Exponentialfunktion
log(v1)         % Nat�rlicher Logarithmus
log10(v1)       % Zehnerlogarithmus

abs(v1)         % Betrag
sign(v1)        % Signum
round(v1)       % Runden
ceil(v1)        % Nach oben runden
floor(v1)       % Nach unten runden

real(2+i)               % Realteil
imag(2+i)               % Imagin�rteil
conj(2+i)               % komplexe Konjugation
angle(2+i)              % Phase einer komplexen Zahl

sin(v1)         % Sinus
cos(v1)         % Cosinus
tan(v1)         % Tangens
cot(v1)         % Cotangens


%%

% Eingabe von Vektoren und Matrizen

zv  = [1,2,3]   % Zeilenvektor 
sv  = [1;2;3]   % Spaltenvektor
zv2 = 1:1:5     % Zeilenvektor [1 2 3 4 5]  
                % Verwendung start:(schrittweite:)ende
                % Der Parameter schrittweite ist optional
zv3 = 1:2:5     % Zeilenvektor [1 3 5]
zv4 = linspace(1,2,11) % Zeilenvektor aus 11 �quidistanten Werten von 1 bis 2

zv5 = [zv zv3]  % Zusammenf�gen von Vektoren

% Eingabe von Matrizen

matr1 = [1 2;3 4] % 2x2 Matrix
matr2 = ones(3)   % 3x3 Einsmatrix
matr3 = zeros(4,2)% 4x2 Nullmatrix
matr4 = randn(4)  % 4x4 Zufallsmatrix
matr5 = rand(4)   % 4x4 Zufallsmatrix mit Werten zwischen 0 und 1
matr6 = eye(3)    % 3x3 Einheitsmatrix 

% Weitere vorgegebene Zahlen, Matrizen und Vektoren sind hier zu finden
% doc elmat              

% Zugriff auf Eintr�ge von Vektoren und Matrizen

zv(2)           % Zweites Element von zv
zv(end)         % Letztes Element von zv
matr1(2,2)      % Element in der zweiten Spalte und der zweiten Zeile
length(zv)      % L�nge des Vektors zv
size(matr1)     % Gr��e der Matrix matr1
matr2(1:2,2:3)  % Auf gewisse Submatrizen zur�ckgreifen

% Rechnen mit Matrizen und Vektoren

det(matr1)      % Determinante der Matrix
rank(matr1)     % Rang der Matrix
eig(matr1)      % Eigenwerte der Matrix
inv(matr1)      % Inverse der Matrix
matr1'          % Transponierte der Matrix

zv*zv'          % Skalarprodukt
zv'*zv          % Dyadisches Produkt
zv.*zv          % Elementweises Ausf�hren der Multiplikation

matr1*[1;2]     % Matrix-Vektor-Produkt
matr1^2         % Matrixquadrat

[min1,min2] = min(zv) % Minimum und Lage des Minimums eines Vektors
[max1,max2] = max(zv) % Maximum und Lage des Maximums eines Vektors
mean(zv)        % Mittelwert eines Vektors
std(zv)         % Standartabweichung vom Mittelwert eines Vektors
sum(zv)         % Summe der Vektorelemente
prod(zv)        % Produkt der Vektorelemente

%%

% Strukturen und Zellen

% Um komplexe zusammengeh�rige Daten �bersichtlich zu verwalten, empfiehlt
% sich die Verwendung von Strukturen. Der Zugriff erfolgt mit Hilfe des
% "."-Operators

% Erzeugen einer Struktur

stru1 = struct('name','Albert','age',21)    % Mit dem Struct-Operator
stru2.name = 'Albert'                       % Mit dem Punkt-Operator
stru2.age = 129

% Der Zugriff auf die Struktur erfolgt analog

% Verschachtelte Struktur

stru3.name.vorname = 'Albert'
stru3.name.nachname = 'Einstein'
stru3.age = 129

% Erzeugen eines Zweiten Eintrages in der Struktur

stru3(2).name.vorname = 'Isaac'
stru3(2).name.nachname = 'Newton'
stru3(2).age = 366


% Cell-Arrays: Werden meist zum Abspeichern unterschiedlicher Daten benutzt

zell1 = cell(2,2)      % Erzeugen der Zelle mit dem cell-Befehl
zell1{1,1} = 1         % Zugriff auf die Elemente mit geschweiften Klammern
zell1{1,2} = [1;2]
zell1{2,1} = 'Hello World'
zell1{2,2} = date

zell2{1,1} = 1         % Direkte Erzeugung ohne vorherige Deklaration
zell2{1,2} = [1;2]
zell2{2,1} = 'Hello World'
zell2{2,2} = date



%%

% Verzweigungen und Schleifen

% If-Verzweigungen

var = 3;

if var > 2
    display('Die Variable ist gr��er als 2') %Ausgabe im Command Window
else
    display('Die Variable ist kleiner als 2 oder gleich 2')
end

% Verschachtelte Verzweigungen

if ((var > 1) && (var<3))
    display('Die Variable liegt zwischen 1 und 3')
elseif (var>=3)
    display('Die Variable ist gr��er als 3 oder gleich 3')
else
    display('Die Variable ist kleiner als 1 oder gleich 1')
end

% Case Anweisung

switch var case 1; a=1, case{3,4,5}; a=5, otherwise a=10, end

% Schleifen

for i=1:2:10
    display(['Der aktuelle Index betr�gt: ',num2str(i)])
    % Umwandlung von Zahlen in Strings mit dem Befehl num2str
    pause(1)     % Einf�hren von Pausen bei der Ausgabe
end

k=0;

while k<5
    display(sprintf('Der aktuelle Index betr�gt: %i',k))
    % Umwandlung von Zahlen in Strings mit dem Befehl num2str
    pause(1)     % Einf�hren von Pausen bei der Ausgabe
    k=k+2;
end

% Das Abbrechen von Schleifen (und m-code allgemein) ist mit dem Befehl break m�glich


while k>0
    display(sprintf('Der aktuelle Index betr�gt: %i',k))
    % Umwandlung von Zahlen in Strings mit dem Befehl num2str
    pause(1)     % Einf�hren von Pausen bei der Ausgabe
    if k>10
        display('Schleife unterbrochen')
        break
    end
    k=k+2;
end

% Anmerkung: Sollte unbeabsichtigt eine Endlosschleife entstanden sein  
%            (oder eine andere Berechnung nicht terminieren), so kann 
%            diese mittels STRG+C unterbrochen werden

%%

% Matlabfunktionen:
% Dies sind spezielle m-files, welchen Parameter �bergeben werden k�nnen
% und dann das Ergebnis der Funktion zur�ckgeben. Der Name der Funktion
% muss immer gleich dem Namen des m-files sein
% Syntax: funtion [out] = name(in)  Abschluss der Funktion mit end

help mittelwert   % Ausgabe der Hilfe

mittelwert(2,3)   % Ausgabe des Ergebnis

% Anmerkung: Variablen in einer Funktion sind lokal, sofern nicht anders
% deklariert. Allerdings sollte darauf wenn m�glich verzichtet werden.

%% 

% Datenexport:
% Daten werden in Matlab standardm��ig als sogenanntes mat-file exportiert

vasa = [1 2 3];

save sav_dat vasa % Abspeichern des Vektors vasa in der Datei sav_dat.mat

vasa2 = [4 5 6];

save sav_dat -append vasa2  %Anh�ngen des Vektors vasa2 an die existierende Datei

clear all

load sav_dat.mat   % Laden einer mat-Datei
vec = load('sav_dat.mat') % Direktes Zuweisen zu einer Variable

% Weitere Datenexporte und Importe

vasa = [1 2 3];

save sav_dat_asc.txt -ascii vasa   %Abspeichern in einer ascii-Datei (z.B. txt)
clear all
load -ascii sav_dat_asc.txt  % Laden der ascii-Datei


vasa = [1 2 3];

xlswrite('sav_dat_xls',vasa);    % Abspeichern in einer Excel-Datei
clear all
var = xlsread('sav_dat_xls')

%%

% Abschluss: Einige zus�tzliche Features

% Erzeugen einer Stoppuhr um die Berechnungszeit auszugeben

tic
toc

tic
pause(1)
toc

% Angabe der Variablen im Arbeitsspeicher

who 
whos  % Ausgabe mit Variablentyp


% L�schen einzelner Variablen

a = 1;
clear a

%% Functions handles

u = @(t) sin(2*t);

t_vec=0:1e-3:4;
u_vec=u(t_vec);
figure();
plot(t_vec,u_vec);


