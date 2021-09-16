%First define all parameters
r = 1; % [day^-1] Growth rate of diatoms - Visser2012
K = 2; % [gC / m3] Carrying capacity of diatoms - Riemann1999
wN = 10^-7; % [gC] Carbon weight of one ciliate
wP = 10^-5; % [gC] Carbon weight of one copepod
bN = 10^6*(3*10^-5)^3; % [m3 / day] Clearance rate of ciliates - 10^6 body volumes per day - Kiorboe 2011
bPR = 10^6*(10^-4)^3; % [m3 / day] Clearance rate of copepods - 10^6 body volumes per day - Kiorboe 2011
bPN = 10^6*(10^-4)^3; % [m3 / day] Clearance rate of copepods - 10^6 body volumes per day - Kiorboe 2011
epsN = 0.8; % [-] Assimilation efficiency for ciliates - Stoecker1984
epsPR = 0.7; % [-] Assimilation efficiency for copepods - Pinti2019b
epsPN = 0.7; 
muR = 0; % [day^-1] Background mortality rate for diatoms
muN = 0.1; % [day^-1] Background mortality rate for ciliates
muP = 0.1; % [day^-1] Background mortality rate for copepods - Visser2007
fcost1 = @(tau,S) (1-S).^tau; % [-] Cost of a defense (in resource acq.) of level S when the tradeoff shape is governed by tau - Cadier2019
fcost2 = @(tau,S) 0;%(1-S).^tau; % [-] Cost of a defense (in extra metabolism) of level S when the tradeoff shape is governed by tau - Cadier2019
fgain = @(tau,S) 1-S.^tau; % [-] Gain for a defense of level S when the tradeoff shape is governed by tau - Cadier2019
tau = 0.5; % [-] Parameter for the trade-off shape - we choose it for now

dRdt = @(C, t, S) r*(1-C(1)/K)*C(1)-fcost1(tau,S)*bN*C(2)*C(1)/wN - bPR*C(1)*C(3)/wP - muR*C(1); % [gC m^-3 day^-1] Variation in diatom concentration over time
dNdt = @(C, t, S) epsN*fcost1(tau,S)*bN*C(2)*C(1)/wN - fcost2(tau,S)*C(2) - fgain(tau,S)*bPN*C(3)*C(2)/wP - muN*C(2);  % [gC m^-3 day^-1] Variation in ciliate concentration over time
dPdt = @(C, t, S) epsPR*bPR*C(1)*C(3)/wP + epsPN*fgain(tau,S)*bPN*C(3)*C(2)/wP - muP*C(3); % [gC m^-3 day^-1] Variation in copepod concentration over time

dCdt = @(C,t,S) [dRdt(C,t,S); dNdt(C,t,S); dPdt(C,t,S)]; % [gC m^-3 day^-1] system of ODEs



%% Solve the system of ODEs
S = .5; % [-] level of defense adopted by the ciliates
options = odeset('NonNegative',[1 2 3]);
[t, C] = ode15s(@(t,C) dCdt(C,t,S), 0:10:1000, [1 .1 .1],options);

figure, 
subplot(311)
plot(t,C(:,1))
ylabel('Diatoms')
hold on

subplot(312)
plot(t,C(:,2))
ylabel('Ciliates')
hold on

subplot(313)
plot(t,C(:,3))
ylabel('Copepods')
hold on

%% Compute detritus creation rate
phi = 0.4; % [-] (1-assim efficiency)/assim efficiency - ratio between growth and what is excreted as FP
R = C(end,1); %  [gC m^-3] equilibrium concentration of diatoms
N = C(end,2); %  [gC m^-3] equilibrium concentration of ciliates
P = C(end,3); %  [gC m^-3] equilibrium concentration of copepods


Crea_R = 0; %phi*r*(1-C(1)/K)*C(1) %assume only leakage from diatoms
Crea_N = (1-epsN)*fcost1(tau,S)*bN*N*R/wN; % [gC m^-3 day^-1] Creation rate of N fecal pellets
Crea_P = (1-epsPR)*bPR*R*P/wP + (1-epsPN)*fgain(tau,S)*bPN*N*P/wP; % [gC m^-3 day^-1] Creation rate of P fecal pellets

Crea_Rc = muR*R; % [gC m^-3 day^-1] Creation rate of carcasses from R
Crea_Nc = muN*N; % [gC m^-3 day^-1] Creation rate of carcasses from N
Crea_Pc = muP*P; % [gC m^-3 day^-1] Creation rate of carcasses from P

%% Now plot the flux
mld = 10; % [m] Mixed layer depth
alpha = @(z) 0.5; % [day^-1] Degradation rate - assumed a Q10 of 2 and a reference temperature of 18

uRfp = 10^-5; % [m/day] Not used, sinking rate of R fecal pellets
uNfp = 0.1; % [m/day] Sinking rate of N fecal pellet
uPfp = 50; % [m/day] Sinking rate of P fecal pellet
uRc = 1; % [m/day] Sinking rate of R carcasses
uNc = 10; % [m/day] Sinking rate of N carcasses
uPc = 500; % [m/day] Sinking rate of P carcasses
speed = [uRfp uNfp uPfp uRc uNc uPc]'; % [m/day] all speeds

Zwc = 0:5:500;
Flux = zeros(size(Zwc,2),1);

So = @(z) (z<mld)*[Crea_R Crea_N Crea_P Crea_Rc Crea_Nc Crea_Pc]';
[z, C] = ode45(@(z,c) (So(z)- alpha(z)*c)./speed,Zwc, [0 0 0 0 0 0]); % [C is gC / m^3 concentration of detritus in the water column]


F = speed'.*C; % [gC / m2 / day] Flux at each depth
Ftot = sum(F,2); % [gC / m2 / day] Total flux at each depth

figure,
plot(Ftot, Zwc)
ylabel('Depth [m]')
xlabel('Flux [gC m^{-2} day^{-1}]')
set(gca,'ydir','reverse')
