
function dndt = ODE_solver(t, substance, initialSubstance, volume, tm)

rateConstants.k1 = 1e11;    %Acid-base reaction
rateConstants.k1r = 5.6e1;  %Backreaction acid-base
rateConstants.k3 = 5.6e9;   % L/mol/s Data from Ruasse et al. 1986
rateConstants.k3r = 7.5e6;  % L/mol/s Data from Ruasse et al. 1986

%linear incorporation function
g = t/tm;
dg = 1/tm;

%volume V2(t)
v = volume.injection + volume.main * g;

%time-dependant ionic strength for perchloric acid and a borate buffer
ionicstrength = 0.5 * (initialSubstance.ClO4 + substance(1) + 2 * substance(2) + 2 * substance(3) + 2 * substance(4) + substance(6)) / v;

%reaction rate constant of the Dushman reaction from Arian and Pauer 2021
rateConstants.k2 = 1.37E9*10^(-1.93*sqrt(ionicstrength)/(1+sqrt(ionicstrength))+0.40*ionicstrength);

%substance(1) = H+; substance(2) = H2BO3-; substance(3) = IO3-;
%substance(4) = I-; substance(5) = I2; substance(6) = I3-; 
%substance(7) = H3BO4;

%reaction rate equation of the acid-base reaction: H2BO3- + H+ -> H3BO3
%for the Dushman reaction:  6H+ + 5I- + IO3- -> 3I2 + 3H2O
%iodine-triiodide equilibrium reaction: I2 + I- -> I3-
r1 = rateConstants.k1 / v * substance(1) * substance(2)-rateConstants.k1r*substance(7);    
r2 = rateConstants.k2 / v^4 * substance(1)^2 * substance(4)^2 *substance(3);    
r3 = rateConstants.k3 / v * substance(4) * substance(5) - rateConstants.k3r * substance(6);

dndt = zeros(length(substance),1);
dndt(1) = -r1 - 6 * r2;                            %dn(H+)/dt
dndt(2) = dg * initialSubstance.H2BO3 - r1;        %dn(H2BO3-)/dt
dndt(3) = dg * initialSubstance.IO3 - r2;          %dn(IO3-)/dt
dndt(4) = dg * initialSubstance.I - 5 * r2 - r3;   %dn(I-)/dt
dndt(5) = 3 * r2 - r3;                             %dn(I2)/dt
dndt(6) = r3;                                      %dn(I3-)/dt
dndt(7) = dg * initialSubstance.H3BO3 + r1;        %dn(H3BO3)/dt

return