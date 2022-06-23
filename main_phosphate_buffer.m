clc
clear all

%volumes V1,0 and V2,0
volume.main =       1.0;	% [L]
volume.injection =	0.001;	% [L]

%For a Phosphate buffer and Perchloric acid
%concentration of each species
concentration.H =       1;                  % concentration H+ in mol/L for strong acid
concentration.ClO4 =    concentration.H;    % concentration ClO4- in case of perchloric acid mol/L
concentration.HPO4=     9.09E-3;            % concentration of HPO4^2- in mol/L
concentration.H2PO4=	9.09E-3;            %concentration of HPO4^- in mol/L
concentration.H3PO4 =   0.0;                %concentration of H3PO4 in mol/L
concentration.I =       5.83E-3;            % concentration of I- in mol/L
concentration.IO3 =     1.17E-3;            % concentration of IO3- in mol/L
concentration.I2 =      0.0;                % concentration of I2 in mol/L
concentration.I3 =      0.0;                % concentration of I3- in mol/L

%initial amount of substances for each species
initialSubstance.H = concentration.H * volume.injection;
initialSubstance.ClO4 = concentration.ClO4 * volume.injection;
initialSubstance.H2PO4 = concentration.H2PO4 * volume.main;
initialSubstance.HPO4 = concentration.HPO4 * volume.main;
initialSubstance.IO3 = concentration.IO3 * volume.main;
initialSubstance.I = concentration.I * volume.main;

%inital amounts of substance
substance = [initialSubstance.H 0 0 0 0 0 0 0];

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

[t,n] = ode15s(@solver_phosphate_buffer, [0 tend], substance, options, initialSubstance, volume, tm);

%Obtained triiodide concentration
I3_end = n(end,8)/(volume.main + volume.injection);




