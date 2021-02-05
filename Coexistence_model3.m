%% Find values of tau and S for which there is coexistence

Npos = @(tau,S) muP*r*bPN/wP*fgain(tau,S)+K*fcost1(tau,S)*bN/wN*bPR/wP*epsN*muP-K*bPR^2/wP^2*epsPR*(muN+fcost2(tau,S))-K*bPR/wP*epsPR*r*bPN/wP*fgain(tau,S); 
Ppos = @(tau,S) K*fcost1(tau,S)*bN/wN*epsN*r*bPN/wP*fgain(tau,S)*epsPN+K*fcost1(tau,S)*bN/wN*bPR/wP*epsPR*(muN+fcost2(tau,S))-K*(fcost1(tau,S)*bN/wN)^2*epsN*muP-(muN+fcost2(tau,S))*r*bPN/wP*fgain(tau,S)*epsPN;

tau = 0:0.1:1;
S = 0:0.1:1;

Npossible = zeros(size(tau,2), size(S,2));
Ppossible = Npossible;

for i=1:size(tau,2) %tau is on the lines
    for j=1:size(S,2) %S is on the columns
        
        Npossible(i,j) = Npos(tau(i),S(j));
        Ppossible(i,j) = Ppos(tau(i),S(j));
        
    end
end

figure,
subplot(121)
imagesc(S,tau,Npossible), colorbar
ylabel('tau')
xlabel('S')
title('Existence of N')
a = max(abs(max(max(Npossible))),abs(min(min(Npossible))));
caxis([-a a])

subplot(122)
imagesc(S,tau,Ppossible), colorbar
xlabel('S')
title('Existence of P')
a = max(abs(max(max(Ppossible))),abs(min(min(Ppossible))));
caxis([-a a])

tot = Npossible > 0 & Ppossible > 0;
figure, imagesc(S,tau,tot), colorbar
xlabel('S')
ylabel('tau')