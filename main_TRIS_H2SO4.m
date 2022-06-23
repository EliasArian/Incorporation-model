clc
clear all

%volumes V1,0 and V2,0
volume.main =       1.0;                % [L]
volume.injection =	0.001;              % [L]

%Calculate initial concentration of HSO4-, SO4^2- and H+ from the
%dissociation of H2SO4
options=odeset('RelTol', 5E-14, 'AbsTol', 5E-16); 
initialConcentration.H2SO4initial  =   0.5;     % concentration H2SO4 in mol/L
Concentration_H2SO4 = [initialConcentration.H2SO4initial 0 0 0];
[t,c] = ode15s(@solver_H2SO4,[0 1],Concentration_H2SO4,options);

%For a TRIS buffer and sulfuric acid
%concentration of each species
concentration.H2SO4  =  c(end,1);               % concentration H2SO4 in mol/L
concentration.HSO4   =  c(end,2);               % concentration HSO4- in mol/L
concentration.SO4 	=   c(end,3);               % concentration SO4^2- in mol/L
concentration.H 	=   c(end,4);               % concentration H+ in mol/L
concentration.TRIS =    9.09E-3;                % concentration TRISH in mol/L
concentration.TRISH =   9.09E-3;                % concentration TRISH+ in mol/L
concentration.ClO4 =    concentration.TRISH;    % concentration ClO4- in mol/L
concentration.I =       5.83E-3;                % concentration of I- in mol/L
concentration.IO3 =     1.17E-3;                % concentration of IO3- in mol/L
concentration.I2 =      0.0;                    % concentration of I2 in mol/L
concentration.I3 =      0.0;                    % concentration of I3- in mol/L

initialSubstance.H2SO4 = concentration.H2SO4 * volume.injection;
initialSubstance.HSO4 = concentration.HSO4 * volume.injection;
initialSubstance.SO4 = concentration.SO4 * volume.injection;
initialSubstance.H = concentration.H * volume.injection;
initialSubstance.TRIS = concentration.TRIS * volume.main;
initialSubstance.TRISH = concentration.TRISH * volume.main;
initialSubstance.ClO4 = concentration.ClO4 * volume.main;
initialSubstance.I = concentration.I * volume.main;
initialSubstance.IO3 = concentration.IO3 * volume.main;

%inital amounts of substance
substance = [initialSubstance.H2SO4 initialSubstance.HSO4 initialSubstance.SO4 initialSubstance.H 0 0 0 0 0 0 0];

%Examined tm
tm = 0.01;

%Interval of integration
%tend is the actual micromixing time which corresponds to te
%tend for Fournier et al.'s linear function; 
%tend = volume.main / volume.injection * tm;

%tend for Fournier et al.'s exponential function; 
tend = log( (volume.main + volume.injection)/ volume.injection ) * tm;

% Setting tolerances for odesolver
options=odeset('RelTol', 5E-14, 'AbsTol', 5E-16); 

[t,n] = ode15s(@solver_TRIS_H2SO4, [0 tend], substance, options, initialSubstance, volume, tm);

%Obtained triiodide concentration
I3_end = n(end,10)/(volume.main+volume.injection);





