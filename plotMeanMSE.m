function plotMeanMSE(meanMSE,BS,user, N)
% This function plots the mean MSE calculated after the Montecarlo method
% meanMSE: matrix of MSE calculated
% BS: The number of the BS
% user: The user in the BS
% N: Vector of Number of antennas at the terminals


figure;
plot(N,10*log10(real(meanMSE(BS,:,user))));
xlabel('N (Antennas at terminals)')
ylabel('MSE(dB)')
title(['MSE of BS ',num2str(BS)]);

end