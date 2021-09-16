%% Plot the trade-off curve
S = linspace(0,1,100);
T = 0.9;
plot(1-fgain(T,S),1-fcost1(T,S),'k')
xlabel('Reduction in predation')
ylabel('Reduction in foraging')
hold on

plot(1-fgain(.5,S),1-fcost1(.5,S),'k')
plot(1-fgain(.2,S),1-fcost1(.2,S),'k')