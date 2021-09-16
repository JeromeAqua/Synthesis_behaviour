phi = 0.4; % [-] (1-assim efficiency)/assim efficiency - ratio between growth and what is excreted as FP

r = 1; % [day^-1] Growth rate
K = 1; % [gC m^-3] Carrying capacity
beta = 8*10^-6; % [m^3 day^-1] Predator clearance rate
wp = 10^-5; % [gC] Predator carbon weight
eps = 0.6; % [-] Predator assimilation efficiency
muP = 0.1; % [day^-1] Predator mortality rate
muN = 0.1; % [day^-1] Prey mortality rate

fcost1 = @(tau,S) (1-S).^tau; % [-] Cost of a defense (in resource acq.) of level S when the tradeoff shape is governed by tau - Cadier2019
fcost2 = @(tau,S) 0;%(1-S).^tau; % [-] Cost of a defense (in extra metabolism) of level S when the tradeoff shape is governed by tau - Cadier2019
fgain = @(tau,S) 1-S.^tau; % [-] Gain for a defense of level S when the tradeoff shape is governed by tau - Cadier2019
tau = 0.5; % [-] Parameter for the trade-off shape - we choose it for now
% S = 0.8;

options = odeset('NonNegative',[1 2]);


S = 0:.01:1;
N= zeros(size(S));
P = zeros(size(S));

for i=1:size(S,2)
    s = S(i);
    dCdt = @(C,t) [r*(1-C(1)/K)*C(1)*fcost1(tau,s) - beta*fgain(tau,s)*C(1)*C(2)/wp - fcost2(tau,s)*C(1) - muN*C(1);...
                   eps*beta*fgain(tau,s)*C(1).*C(2)/wp - muP*C(2)            ];% [gC m^-3 day^-1] change of pop per time
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
subplot(221)
plot(S,N)
ylabel({'Concentration','[gC m^{-3}]'})
title('prey')
ylim([0 1.2])
subplot(222)
plot(S,P)
title('predator')
% ylabel('Predator [gC m^{-3}]')
ylim([0 1.2])

    

Crea_N = phi*r*(1-N/K).*N.*fcost1(tau,S); % [gC m^-3 day^-1] Creation rate of detritus from prey fecal pellets
Crea_P = (1-eps)*beta*fgain(tau,S).*N.*P/wp; % [gC m^-3 day^-1] Creation rate of detritus from predator fecal pellets

Crea_Nc = muN*N; % [gC m^-3 day^-1] Creation rate of carcasses from prey
Crea_Pc = muP*P; % [gC m^-3 day^-1] Creation rate of carcasses from predators

subplot(223)
plot(S,Crea_N)
ylabel({'Detritus production','rate [gC m^-^3 day^-^1]'})
hold on
plot(S,Crea_Nc)
plot(S,Crea_Nc+Crea_N,'k')
ylim([0 .18])

subplot(224)
plot(S,Crea_P)
hold on
plot(S,Crea_Pc)
plot(S,Crea_Pc+Crea_P,'k')
ylim([0 .18])


%% Now plot the flux
mld = 10; % [m] Mixed layer depth
T = @(z) 4 + 20*(1 - tanh((-mld +z)/100))/2; % [deg C] Temperature
alpha = @(z) 0.5; %*2.^((T(z)-18)/10); % [day^-1] Degradation rate - assumed a Q10 of 2 and a reference temperature of 18

uNfp = @(s) linspace(0.1,0.1,size(s,2)); % [m/day] Sinking rate of prey fecal pellet
uPfp = @(s) (s+1)*50; % [m/day] Sinking rate of prey fecal pellet
uNc = @(s) (s+1)*10; % [m/day] Sinking rate of prey carcasses
uPc = @(s) linspace(500,500,size(s,2)); % [m/day] Sinking rate of predator carcasses

SPEED = @(s) [uNfp(s); uPfp(s); uNc(s); uPc(s)]; % [m/day] all speeds

Zwc = 0:10:500;
Flux = zeros(size(Zwc,2),size(S,2));
FluxC = zeros(size(Zwc,2),size(S,2));
FluxP = zeros(size(Zwc,2),size(S,2));

for i = 1:size(S,2)
    speed = SPEED(S(i));
    So = @(z) (z<mld)*[Crea_N(i) Crea_P(i) Crea_Nc(i) Crea_Pc(i)]';
    [z, C] = ode45(@(z,c) (So(z)- alpha(z)*c)./speed,Zwc, [0 0 0 0]); % [C is gC / m^3 concentration of detritus in the water column]
    
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

% figure, imagesc(S,Zwc,Flux)
% set(gca,'ydir','reverse')
% hold on
% colorbar
% [~,sopt] = max(Flux'); Sopt = S(sopt); % [-] Optimal defense level for flux at each depth
% plot(Sopt,Zwc,'r')

% figure, 
% plot(S, FluxC(6,:)./Flux(6,:))
% hold on
% plot(S, FluxP(6,:)./Flux(6,:))
% legend('prop prey', 'prop pred')

figure
subplot(132)
plot(S, FluxC(6,:))
hold on
plot(S,FluxP(6,:))
plot(S,Flux(6,:))

subplot(131)
plot(S, FluxC(3,:))
hold on
plot(S,FluxP(3,:))
plot(S,Flux(3,:))

subplot(133)
plot(S, FluxC(11,:))
hold on
plot(S,FluxP(11,:))
plot(S,Flux(11,:))
