function dcdt = solver_H2SO4 (t, Concentration_H2SO4)

ionstrength = 0.5 * (Concentration_H2SO4(4) +Concentration_H2SO4(2) + 4*Concentration_H2SO4(3));

equilibriumConstants.K1 = 1E3;      %EquilibriumConstant H2SO4 -> H+ + HSO4-
% equilibriumConstants.K2 = 10^(log10(0.01023)); %pKs = 1.99
% equilibriumConstants.K2 = 10^(log10(0.01023) + 2.046*sqrt(ionstrength)); %Debye
% equilibriumConstants.K2 = 10^(log10(0.01023) + (2.046*sqrt(ionstrength)/(1+1.5*sqrt(ionstrength))));  %Debye Hückel EquilibriumConstant HSO4- -> H+ + SO4^2-
% equilibriumConstants.K2 = 10^(log10(0.01020) + (2.036*sqrt(ionstrength)/(1+0.719*sqrt(ionstrength)))); %Baes
equilibriumConstants.K2 = 10^(log10(0.01028) + (2.032*sqrt(ionstrength)/(1+0.94*sqrt(ionstrength)))); %Marshal
% equilibriumConstants.K2 = 10^(log10(0.01028) + (2.032*sqrt(ionstrength)/(1+0.512*sqrt(ionstrength)))); %Clegg
% equilibriumConstants.K2 = 10^(log10(0.01050) + (2.036*sqrt(ionstrength)/(1+0.630*sqrt(ionstrength))) + 0.0866 * ionstrength); %Reviewer1


rateConstants.k1 = 1e11;            %Rate constant H+ + HSO4- -> H2SO4
rateConstants.k1r = rateConstants.k1*equilibriumConstants.K1; %Rate constant H2SO4 -> H+ + HSO4-

rateConstants.k2 = 1e11;            %Rate constant H+ + SO4^2- -> HSO4-
rateConstants.k2r = rateConstants.k2*equilibriumConstants.K2; %Rate constant  HSO4- -> H+ + SO4^2-

%H2SO4 -> H+ + HSO4-
r1 = rateConstants.k1r * Concentration_H2SO4(1) - rateConstants.k1 * Concentration_H2SO4(2) * Concentration_H2SO4(4);
%HSO4- -> H+ + SO4^2-
r2 = rateConstants.k2r * Concentration_H2SO4(2) - rateConstants.k2 * Concentration_H2SO4(3) * Concentration_H2SO4(4);

dcdt = zeros(length(Concentration_H2SO4),1);
dcdt(1) = -r1;                      %dH2SO4/dt
dcdt(2) = r1 - r2;                  %dHSO4/dt
dcdt(3) = r2;                       %dSO4/dt
dcdt(4) = r1 + r2;                  %dH/dt

return