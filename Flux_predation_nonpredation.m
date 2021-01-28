% Plot to compute the vertical flux of a phytoplankton dying of natural
% cause of sinking as fecal pellets because of predation

mld = 50; % [m] Mixed layer depth

T = @(z) 4 + 20*(1 - tanh((-mld +z)/100))/2; % [deg C] Temperature
alpha = @(z) 0.5*2.^((T(z)-18)/10); % [day^-1] Degradation rate - assumed a Q10 of 2 and a reference temperature of 18


w1 = 10; % [m/day] Sinking rate of a dead phytoplankton cell
w2 = 100; % [m/day] Sinking rate of a zooplankton fecal pellet

epsZ = 0.7; % [-] Assimilation efficiency of zooplankton

S = @(z) (1 - tanh((-mld +z)/20))/2; % [gC / m3 / day] Source term

[z1, C1] = ode45(@(z,c) (S(z)- alpha(z).*c)./w1,0:10:500, 0);
[z2, C2] = ode45(@(z,c) (S(z)*(1-epsZ)- alpha(z).*c)./w2,0:10:500, 0);

F1 = w1.*C1; % @(z) C1(z).*w1; % [gC / m2 / day] Sinking flux of dead phytoplankton
F2 = w2.*C2; %@(z) C2(z).*w2; % [gC / m2 / day] Sinking flux of fecal pellets

Z = 0:10:500;

figure
subplot(131)
plot(T(Z),Z,'k')
xlabel('Temperature [°C]')
ylabel('Depth [m]')
set(gca,'ydir','reverse')
xlim([0 max(T(Z))])

subplot(132)
plot(S(Z),Z,'k')
xlabel('Source term [gC / m^3 / day]')
set(gca,'ydir','reverse')


subplot(133)
plot(F1,z1,'b')
hold on
plot(F2,z2,'r')
legend('non predation', 'predation')
xlabel('Flux [gC / m^2 / day]')
set(gca,'ydir','reverse')

