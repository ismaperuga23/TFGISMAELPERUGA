function plotMSE(MSE,BS,user,p)
% This function plots the MSE
% MSE: matrix of MSE calculated
% BS: The number of the BS
% user: The user in the BS
% p: vector of power of the pilots

figure;
plot(10*log10(p),10*log10(real(MSE(BS,:,user))));
xlabel('SNR (dB)')
ylabel('MSE(dB)')


end