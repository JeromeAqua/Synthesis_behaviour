phi = 0.4; % [-] (1-assim efficiency)/assim efficiency - ratio between growth and what is excreted as FP
r = 1; % [day^-1] Growth rate
K = 5; % [gC m^-3] Carrying capacity
beta = 8*10^-6; % [m^3 day^-1] Predator clearance rate
wp = 10^-5; % [gC] Predator carbon weight
eps = 0.6; % [-] Predator assimilation efficiency
m = 0.1; % [day^-1] Predator mortality rate
mu0 = 0.1; % [day^-1] Prey mortality rate

cost = @(d) (exp(d)-1)./(exp(1)-1); % [day^-1] Defense cost as a function of defense investment
options = odeset('NonNegative',[1 2]);


D = 0:.1:1;
N= zeros(size(D));
P = zeros(size(D));

for i=1:size(D,2)
    d = D(i);
    dCdt = @(C,t) [r*(1-C(1)/K)*C(1) - beta*(1-d)*C(1)*C(2)/wp - cost(d)*C(1) - mu0*C(1);...
                   eps*beta*(1-d)*C(1).*C(2)/wp - m*C(2)            ];% [gC m^-3 day^-1] change of pop per time
    [t, C1] = ode15s(@(t,C) dCdt(C,t), 0:10:10000, [.1 10^-5],options);
% figure, subplot(311)
% plot(t,C1(:,1))
% ylabel('Prey [gC m^{-3}]')
% subplot(312)
% plot(t,C1(:,2))
% ylabel('Predator [gC m^{-3}]')

    N(i) = mean(C1(end-100:end,1));
    P(i) = mean(C1(end-100:end,2));
end

figure, 
subplot(321)
plot(D,N)
ylabel({'Concentration','[gC m^{-3}]'})
title('prey')
subplot(322)
plot(D,P)
title('predator')
% ylabel('Predator [gC m^{-3}]')

    

Crea_N = phi*r*(1-N/K).*N; % [gC m^-3 day^-1] Creation rate of detritus from prey fecal pellets
Crea_P = (1-eps)*beta*(1-D).*N.*P/wp; % [gC m^-3 day^-1] Creation rate of detritus from predator fecal pellets

Crea_Nc = mu0*N; % [gC m^-3 day^-1] Creation rate of carcasses from prey
Crea_Pc = m*P; % [gC m^-3 day^-1] Creation rate of carcasses from predators

subplot(323)
plot(D,Crea_N)
ylabel({'FP source term','[gC m^-^3 day^-^1]'})
subplot(324)
plot(D,Crea_P)

subplot(325)
plot(D,Crea_Nc)
ylabel({'Carcasses source term','[gC m^-^3 day^-^1]'})
subplot(326)
plot(D,Crea_Pc)

%% Now plot the flux
mld = 50; % [m] Mixed layer depth
T = @(z) 4 + 20*(1 - tanh((-mld +z)/100))/2; % [deg C] Temperature
alpha = @(z) 0.5*2.^((T(z)-18)/10); % [day^-1] Degradation rate - assumed a Q10 of 2 and a reference temperature of 18

uNfp = 0.1; % [m/day] Sinking rate of prey fecal pellet
uPfp = 100; % [m/day] Sinking rate of prey fecal pellet
uNc = 10; % [m/day] Sinking rate of prey carcasses
uPc = 500; % [m/day] Sinking rate of predator carcasses

speed = [uNfp uPfp uNc uPc]'; % [m/day] all speeds

Zwc = 0:10:500;
Flux = zeros(size(Zwc,2),size(D,2));
FluxC = zeros(size(Zwc,2),size(D,2));
FluxP = zeros(size(Zwc,2),size(D,2));

for i = 1:size(D,2)
    S = @(z) (z<50)*[Crea_N(i) Crea_P(i) Crea_Nc(i) Crea_Pc(i)]';
    [z, C] = ode45(@(z,c) (S(z)- alpha(z)*c)./speed,Zwc, [0 0 0 0]); % [C is gC / m^3 concentration of detritus in the water column]
    
    F = speed'.*C; % [gC / m2 / day] Flux at each depth
    Ftot = sum(F,2); % [gC / m2 / day] Total flux at each depth
    
    Flux(:,i) = Ftot; % [gC / m2 / day] Total flux at each depth
    FluxC(:,i) = sum(F(:,[1 3]),2); % [gC / m2 / day] 
    FluxP(:,i) = sum(F(:,[2 4]),2); % [gC / m2 / day] 
    
    
%     S = @(z) (z<50)*Crea_P(i);
%     [z2, CPfp] = ode45(@(z,c) (S(z)- alpha(z).*c)./uPfp,0:10:500, 0);
%     
%     S = @(z) (z<50)*Crea_Nc(i);
%     [z3, CNc] = ode45(@(z,c) (S(z)- alpha(z).*c)./uNc,0:10:500, 0);
%     
%     S = @(z) (z<50)*Crea_Pc(i);
%     [z4, CPc] = ode45(@(z,c) (S(z)- alpha(z).*c)./uPc,0:10:500, 0);

end

figure, imagesc(D,Zwc,Flux)
set(gca,'ydir','reverse')
hold on

[~,dopt] = max(Flux'); Dopt = D(dopt); % [-] Optimal defense level for flux at each depth
plot(Dopt,Zwc,'r')

figure, 
plot(D, FluxC(6,:)./Flux(6,:))
hold on
plot(D, FluxP(6,:)./Flux(6,:))
legend('prop prey', 'prop pred')
