function GSM_S(block)

% Simulationsmodell für das elektrische Netzwerk
%
% -------------------------------------------------------------------------
%
% Beschreibung: Simulationsmodell der Gleichstrommaschine mit Propeller
%
% -------------------------------------------------------------------------
%
% Eingaenge:    u1(1) ... uGSM   Eingangsspannung [V]
%               u1(2) ... Mext   externes Moment [Nm]
%
% Zustaende:    x(1)  ... iGSM
%               x(2)  ... phiGSMP
%               x(3)  ... wGSM
%               x(4)  ... wP
%
% Ausgaenge:    y1(1) ... iGSM
%               y1(2) ... phiGSMP
%               y1(3) ... wGSM
%               y1(4) ... wP
%
%               y2(1) ... MGSM
%               y2(1) ... Mkopp
%
% Parameter: 
% iGSM_0, phiGSMP_0, wGSM_0, wP_0, LGSM, RGSM, kGSM, JGSM, dcGSM, dvGSM, JP, dcP, dvP, dqP, cGSMP, dGSMP
%
%               p(1) ... iGSM_0
%               p(2) ... phiGSMP_0
%               p(3) ... wGSM_0
%               p(4) ... wP_0
%               p(5) ... LGSM = 1.4e-3;
%               p(6) ... RGSM = 0.46;
%               p(7) ... kGSM = 0.1;
%               p(8) ... JGSM = 12.4e-3;
%               p(9) ... dcGSM = 0.152;
%               p(10) ... dvGSM = 1.8e-3;
%               p(11) ... JP = 32.5e-3;
%               p(12) ... dcP = 0.169;
%               p(13) ... dvP = 2.7e-3;
%               p(14) ... dqP = 1e-4;
%               p(15) ... cGSMP = 0.6822;
%               p(16) ... dGSMP = 1e-5;
%
% -------------------------------------------------------------------------
% Abtastzeit (sample time): zeitkontinuierlich (continuous)
% -------------------------------------------------------------------------


% Die Funktion setup (s.u.) dient der Initialiserung des Matlab Objektes
% (block). Im Objekt block sind alle fuer die Simulation in Simulink
% notwendigen Eigenschaften (Eingaenge, Zustaende, Ausgaenge, Parameter,
% usw.) des dynamischen Systems (math. Modell) zusammengefasst.
setup(block);

% -------------------------------------------------------------------------
% Initialisierung des Simulationsobjektes block
% -------------------------------------------------------------------------

function setup(block)
  
  % Anzahl der Eingangs- und Ausgangsports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 2;
  
  % Anzahl der zeitkontinuierlichen Zustaende
  block.NumContStates = 4;

  % Anzahl der Parameter
  block.NumDialogPrms = 16;
  
  % Dimensionen der Eingangsports
  % Flag DirectFeedthrough kennzeichnet, ob ein Eingang direkt an einem
  % Ausgang auftritt, d.h. y=f(u)
  block.InputPort(1).Dimensions        = 2;
  block.InputPort(1).SamplingMode = 'Sample';
  block.InputPort(1).DirectFeedthrough = false;

  % Dimensionen der Ausgangsports  
  block.OutputPort(1).Dimensions       = 4;
  block.OutputPort(1).SamplingMode = 'Sample';
  block.OutputPort(2).Dimensions       = 2;
  block.OutputPort(2).SamplingMode = 'Sample';
  
  % Einstellen der Abtastzeit: [0 0] wird verwendet fuer die
  % zeitkontinuierliche Simulation.
  block.SampleTimes = [0 0];
  
  % ------------------------------------------------
  % NICHT VERAENDERN
  % ------------------------------------------------
  % 
  % Registrieren der einzelnen Methoden
  % Hier: InitializeConditions ... Initialisierung
  %       Outputs ...       Berechnung der Ausgaenge
  %       Derivatives ...   Berechnung der Zustaende
  %       Terminate ...     Konsistentes Beenden der Simulation

  block.RegBlockMethod('InitializeConditions',    @InitConditions); 
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivatives);  
  block.RegBlockMethod('Terminate',               @Terminate);


% -------------------------------------------------------------------------
% Setzen der Anfangsbedingungen der Zustaende
% -------------------------------------------------------------------------

function InitConditions(block)
  
  % Einlesen der Parameter des Systems
  iGSM_0  = block.DialogPrm(1).Data;
  phiGSMP_0 = block.DialogPrm(2).Data;
  wGSM_0  = block.DialogPrm(3).Data;
  wP_0    = block.DialogPrm(4).Data;
  
  % Eingabe der Anfangsbedingungen
  x0(1) = iGSM_0;
  x0(2) = phiGSMP_0;
  x0(3) = wGSM_0;
  x0(4) = wP_0;
  
  % Schreiben auf Objekt block (NICHT VERAENDERN)
  block.ContStates.Data = x0;


% -------------------------------------------------------------------------
% Berechnen der Ausgaenge
% -------------------------------------------------------------------------

function Output(block)
   
  % Einlesen der Parameter des Systems
  kGSM  = block.DialogPrm(7).Data;
  cGSMP = block.DialogPrm(15).Data;
  dGSMP = block.DialogPrm(16).Data;

  % Shortcut fuer die Zustaende
  x = block.ContStates.Data;
  iGSM    = x(1);
  phiGSMP = x(2);
  wGSM    = x(3);
  wP      = x(4);
  
  % Berechnung der Ausgaenge
  % Port 1:  
  y1(1) = iGSM;
  y1(2) = phiGSMP;
  y1(3) = wGSM;
	y1(4) = wP;
  
  % Berechnung der Ausgaenge
  % Port 2:
  Mkopp = (wGSM-wP)*dGSMP + (phiGSMP)*cGSMP;
  MGSM = kGSM*iGSM;
  
  y2(1) = MGSM;
  y2(2) = Mkopp;
  
  % Schreiben auf Objekt block
  block.OutputPort(1).Data = y1;
  block.OutputPort(2).Data = y2;

% -------------------------------------------------------------------------
% Berechnen der Zustaende
% -------------------------------------------------------------------------

function Derivatives(block)

  % Einlesen der Parameter des Systems
  LGSM  = block.DialogPrm(5).Data;
  RGSM  = block.DialogPrm(6).Data;
  kGSM  = block.DialogPrm(7).Data;
  JGSM  = block.DialogPrm(8).Data;
  dcGSM = block.DialogPrm(9).Data;
  dvGSM = block.DialogPrm(10).Data;
  JP    = block.DialogPrm(11).Data;
  dcP   = block.DialogPrm(12).Data;
  dvP   = block.DialogPrm(13).Data;
  dqP   = block.DialogPrm(14).Data;
  cGSMP = block.DialogPrm(15).Data;
  dGSMP = block.DialogPrm(16).Data;
  
  % Shortcut fuer den Eingang
  u = block.InputPort(1).Data;
  uGSM = u(1);
  Mext = u(2);
  
  % Shortcut fuer die Zustaende
  x = block.ContStates.Data;
  iGSM    = x(1);
  phiGSMP = x(2);
  wGSM    = x(3);
  wP      = x(4);
  
  %Gleichungen
  MP = dcP + dvP*wP+dqP*wP^2+Mext;
  MrGSM = dcGSM + dvGSM*wGSM;
  Mkopp = (wGSM-wP)*dGSMP + (phiGSMP)*cGSMP;
  MGSM = kGSM*iGSM;
  
  % Berechnen der Zeitableitungen der Zustaende
  dx(1) = (1/LGSM)*(uGSM-RGSM*iGSM-kGSM*wGSM);
  dx(2) = wGSM - wP;
  dx(3) = (MGSM-MrGSM-Mkopp)/(JGSM);
  dx(4) = (Mkopp-MP)/JP;
  
  % Schreiben auf Objekt block
  block.Derivatives.Data = dx;


% -------------------------------------------------------------------------
% Operationen am Ende der Simulation
% -------------------------------------------------------------------------

% Die function Terminate wird hier nicht verwendet,
% muss aber vorhanden sein!
function Terminate(block)

