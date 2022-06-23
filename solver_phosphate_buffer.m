
function dndt = solver_phosphate_buffer(t, substance, initialSubstance, volume, tm)

%exponential incorporation function according to Fournier et al.
g = exp(t/tm);
dgdt = exp(t/tm)/tm;

%linear incorporation function according to Fournier et al.
% g = 1+t/tm;
% dgdt = 1/tm;

%volume growth of V2(t)
v = volume.injection * g;

%time-dependant ionic strength for perchloric acid and a phosphate buffer
ionstrength = 0.5 * (initialSubstance.ClO4 + substance(1) + 6*substance(2) + 2*substance(3) + 2*substance(4) + substance(6)+ 2*substance(7)) / v;

%ionic strength dependency of pK2 for phosphoric acid
%equilibriumConstants.pK2 = 7.037-1.51*sqrt(ionstrength)/(1+0.329*10^8*5.11*10^-8*sqrt(ionstrength));   %Shuma
equilibriumConstants.pK2 = -log10(6.31E-8)- 1.53*sqrt(ionstrength)/(1+1.48*sqrt(ionstrength));          %extended Debye
%equilibriumConstants.pK2 = 7.16 - 1.5*sqrt(ionstrength)/(1+1.65*sqrt(ionstrength));                    %Cohn
%equilibriumConstants.pK2 = 7.2;                                                                        %fixed pKa
equilibriumConstants.K2 = 10^(-equilibriumConstants.pK2);   %H2PO4- -> HPO4^2- + H+
equilibriumConstants.K1 = 7.4E-3;                           %H3PO4 -> H2PO4- + H+

rateConstants.k1 = 1E11;                                        %HPO4^2- + H+ -> H2PO4-
rateConstants.k1r = rateConstants.k1*equilibriumConstants.K2;   %H2PO4- -> HPO4^2- + H+
rateConstants.k2 = 1.32E11;                                     %H2PO4- + H+ -> H3PO4 
rateConstants.k2r = 3.8E8;                                      %H3PO4 + H+ -> H2PO4-

rateConstants.k4 = 5.6e9;   % L/mol/s Data from Ruasse et al. 1986
rateConstants.k4r = 7.5e6;  % 1/s Data from Ruasse et al. 1986

%reaction rate constant of the Dushman reaction from Arian and Pauer 2021
rateConstants.k3 = 1.37E9*10^(-1.93*sqrt(ionstrength)/(1+sqrt(ionstrength))+0.40*ionstrength);

%substance(1) = H+; substance(2) = HPO4^2-; substance(3) = H2PO4^-;
%substance(4) = H3PO4; substance(5) = IO3-; substance(6) = I-;  
%substance(7) = I2; substance(8) = I3; substance(9) = V2;

%reaction rate equation of the acid-base reaction: HPO4^2- + H+ -> H2PO4^-
%reaction rate equation of the acid-base reaction: H2PO4^- + H+ -> H3PO4
%for the Dushman reaction:  6H+ + 5I- + IO3- -> 3I2 + 3H2O
%iodine-triiodide equilibrium reaction: I2 + I- -> I3
r1 = rateConstants.k1 / v * substance(1) * substance(2) - rateConstants.k1r * substance(3);
r2 = rateConstants.k2 / v * substance(1) * substance(3) - rateConstants.k2r * substance(4);
r3 = rateConstants.k3 / v^4 * substance(1)^2 * substance(6)^2 * substance(5);
r4 = rateConstants.k4 /v * substance(6) * substance(7) - rateConstants.k4r * substance(8);

dndt = zeros(length(substance),1);
dndt(1) = -r1 - r2 - 6 * r3;                                                        %dn(H+)/dt
dndt(2) = dgdt * volume.injection/volume.main * initialSubstance.HPO4 - r1;         %dn(HPO4^2-)/dt
dndt(3) = dgdt * volume.injection/volume.main  * initialSubstance.H2PO4 + r1 - r2;  %dn(H2PO4^-)/dt
dndt(4) = r2;                                                                       %dn(H3PO4)/dt
dndt(5) = dgdt * volume.injection/volume.main  * initialSubstance.IO3 - r3;         %dn(IO3-)/dt
dndt(6) = dgdt * volume.injection/volume.main  * initialSubstance.I - 5 * r3 - r4;  %dn(I-)/dt
dndt(7) = 3 * r3 - r4;                                                              %dn(I2)/dt
dndt(8) = r4;                                                                       %dn(I3-)/dt
return