clc
clear all

%volumes V1,0 and V2,0
volume.main =       1.0;	% volume bulk V1,0 in L
volume.injection =	0.004;	% volume acid V2,0 in L

%concentration of each species
concentration.H =       1;                  % concentration H+ in mol/L for strong acid
concentration.ClO4 =    concentration.H;    % concentration ClO4- in case of perchloric acid mol/L
concentration.H3BO3 =   0.0909;             % concentration of H3BO3 in mol/L
concentration.H2BO3 =   0.0909;             % concentration of H2BO3- in mol/L
concentration.IO3 =     2.33E-3;            % concentration of IO3- in mol/L
concentration.I =       0.0117;             % concentration of I- in mol/L
concentration.I2 =      0.0;                % concentration of I2 in mol/L
concentration.I3 =      0.0;                % concentration of I3- in mol/L

%initial amount of substances for each species
initialSubstance.H = concentration.H * volume.injection;    %in mol
initialSubstance.ClO4 = initialSubstance.H;                 %in mol
initialSubstance.H3BO3 = concentration.H3BO3 * volume.main; %in mol
initialSubstance.H2BO3 = concentration.H2BO3 * volume.main; %in mol
initialSubstance.IO3 = concentration.IO3 * volume.main;     %in mol
initialSubstance.I = concentration.I * volume.main;         %in mol

%maximum yield of the Dushman reaction (equivalent to maximum segregation)
Yst = 6 * initialSubstance.IO3 / (6 * initialSubstance.IO3 + initialSubstance.H2BO3);

%Setting tolerances for ODE solver;
options=odeset('RelTol', 1E-20, 'AbsTol', 1E-21); 

%Initial conditions
substance = [initialSubstance.H 0 0 0 0 0 0];

%Examined micromixing time tm
tm = 0.01;

%Interval of integration
%linear function
tend = tm;
%exponential function at 5*tm
%tend = 5*tm;

%ODE solver
[t,n] = ode15s(@ODE_solver, [0 tend], substance, options, initialSubstance, volume, tm);

%Obtained Xs from examined tm
Xs = 2*(n(end,5) + n(end,6)) / (initialSubstance.H) / Yst;




