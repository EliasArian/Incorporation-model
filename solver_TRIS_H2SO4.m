
function dndt = solver_TRIS_H2SO4(t, substance, initialSubstance, volume, tm)

g = exp(t/tm);
dgdt = exp(t/tm)/tm;

%Fournier
v = volume.injection * g;

%for sulfuric acid and a TRIS buffer
ionstrength = 0.5 * (substance(2) + 4 * substance(3) + substance(4) + 2*substance(6) + 2 * substance(7) + 2 * substance(8) + substance(10)) / v;

equilibriumConstants.K1 = 1E3;      %EquilibriumConstant H2SO4 -> H+ + HSO4-
% equilibriumConstants.K2 = 10^(log10(0.01023)); %pKs = 1.99
% equilibriumConstants.K2 = 10^(log10(0.01023) + 2.046*sqrt(ionstrength)); %Debye
% equilibriumConstants.K2 = 10^(log10(0.01023) + (2.046*sqrt(ionstrength)/(1+1.5*sqrt(ionstrength))));  %Debye Hückel
% equilibriumConstants.K2 = 10^(log10(0.01020) + (2.036*sqrt(ionstrength)/(1+0.719*sqrt(ionstrength)))); %Baes
% equilibriumConstants.K2 = 10^(log10(0.01028) + (2.032*sqrt(ionstrength)/(1+0.94*sqrt(ionstrength)))); %Marshal
equilibriumConstants.K2 = 10^(log10(0.01050) + (2.036*sqrt(ionstrength)/(1+0.630*sqrt(ionstrength))) + 0.0866 * ionstrength); %Clegg
equilibriumConstants.K3 = 8.06;     %pKs TRIS

rateConstants.k1r = 1e11;                                       %H+ + HSO4- -> H2SO4
rateConstants.k1 = rateConstants.k1r*equilibriumConstants.K1;   %H2SO4 -> H+ + HSO4-

rateConstants.k2r = 1e11;                                       %H+ + SO4^2- -> HSO4-
rateConstants.k2 = rateConstants.k2r*equilibriumConstants.K2;   %HSO4- -> H+ + SO4^2-

rateConstants.k3 = 1e11;                                            %H+ + TRIS -> TRISH+
rateConstants.k3r = rateConstants.k3*10^(-equilibriumConstants.K3); %TRISH+ -> H+ + TRIS

rateConstants.k5 = 5.6e9;   %L/mol/s Data from Ruasse et al. 1986
rateConstants.k5r = 7.5e6;  %1/s Data from Ruasse et al. 1986

%reaction rate constant of the Dushman reaction from Arian and Pauer 2021
rateConstants.k4 = 1.37E9*10^(-1.93*sqrt(ionstrength)/(1+sqrt(ionstrength))+0.40*ionstrength);

%substance(1) = H2SO4; substance(2) = HSO4-; substance(3) = SO4^2-;
%substance(4) = H+; substance(5) = TRIS; substance(6) = TRISH;  
%substance(7) = I; substance(8) = IO3; substance(9) = I2; substance(10) = I3;

%reaction rate equation for the dissociation of H2SO4: H2SO4 -> H+ + HSO4-
%reaction rate equation for the dissociation of HSO4-: HSO4- -> H+ + SO4^2-
%reaction rate equation for the acid-base reaction: H+ + TRIS -> TRIS-H+
%for the Dushman reaction:  6H+ + 5I- + IO3- -> 3I2 + 3H2O
%iodine-triiodide equilibrium reaction: I2 + I- -> I3


r1 = rateConstants.k1 * substance(1) - rateConstants.k1r / v * substance (2) * substance(4);
r2 = rateConstants.k2 * substance(2) - rateConstants.k2r / v * substance (3) * substance(4);
r3 = rateConstants.k3 / v * substance(5)*substance(4) - rateConstants.k3r * substance(6); 
r4 = rateConstants.k4 / v^4 * substance(4)^2 * substance(7)^2 *substance(8);                 
r5 = rateConstants.k5 / v * substance(7) * substance(9) - rateConstants.k5r * substance(10);

dndt = zeros(length(substance),1);
dndt(1) = -r1;                                                                      %dn(H2SO4)/dt
dndt(2) = r1 - r2;                                                                  %dn(HSO4-)/dt
dndt(3) = r2;                                                                       %dn(SO4^2-)/dt
dndt(4) = r1 + r2 - r3 -  6*r4;                                                     %dn(H+)/dt
dndt(5) = dgdt * volume.injection / volume.main * initialSubstance.TRIS - r3;       %dn(TRIS)/dt
dndt(6) = dgdt * volume.injection / volume.main * initialSubstance.TRISH + r3;      %dn(TRISH+)/dt
dndt(7) = dgdt * volume.injection / volume.main * initialSubstance.I - 5 * r4 - r5; %dn(I-)/dt
dndt(8) = dgdt * volume.injection / volume.main * initialSubstance.IO3 - r4;        %dn(IO3-)/dt
dndt(9) = 3 * r4 - r5;                                                              %dn(I2)/dt
dndt(10) = r5;                                                                      %dn(I3)/dt
dndt(11) = dgdt * volume.injection / volume.main * initialSubstance.ClO4;           %dn(ClO4-)/dt
return