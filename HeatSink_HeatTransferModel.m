% Taken From https://www.sciencedirect.com/science/article/pii/S1359431122013096
% Assumptions: heat transfer to ambient environment is negligible, physical
% properties of the material(density, sp heat, thermal conductivity) are
% not affected by temperature.

% Chamber Parameters
At = 1; % throat area
Dt = (sqrt(At/pi)); %throat diameter 
A = 1; % local cross sectional area
Pc = 1; % chamber pressure
cstar = 1;
R = 1; % throat radius of curvature
T0 = 1;
Tst = 1; % stagnation temp
mu = 1; % viscosity
delta = 1; % wall thickness [m]

cp = 1; % sp. heat capacity 
lamda = 1; % thermal conduct coefficent of wall
dens = 1; %density of wall
Pr = 1; % Prandtl number
mu1 = 1; % eigenvalue
C = 0.026; % heat flow correction coefficent
Twg = 1;

sigma = ((Twg/(2*Tst))*(1+((gamma-1)/2*M*a^2))+0.5)^-0.68 * (1+(gamma-1)/2 *M*a^2)^-0.12;
hg = (C / Dt^0.2)*(mu^0.2*cp/Pr)*(Pc^0.8/cstar)*(Dt^0.1/R)*(At^0.9/A)*sigma;

dT = lamda/(dens*cp)*Twg; % 2.2 - needs to be reviewed/redone

Bi = hg*delta/lamda;

A1 = 2*(sin(mu1)/(mu1+sin(mu1)*cos(mu1)));
B1 = sin(mu1)/mu1;
q = (1-A1*(exp(-mu1^2*F0))*B1*dens*cp*delta*(T0-Tst));
