function plotMSE(MSE,BS,user)
% This function plots the MSE
% MSE: matrix of MSE calculated
% BS: The number of the BS
% user: The user in the BS

figure;
plot(abs(MSE(BS,:,user)));


end