%% In this script we are comparing the effects of the different beamforming
% techniques with 2 UE, letting one be in a fixed position (the one that is
% served by the BS) and the other rotating its angle of arrival.
% UE1 is the one that is served by the BS. UE2 is the one creating the
% interferences

clear all

Nrealizations = 1;
SNR = 5; % value of the fixed SNR in dB (power of noise = 1)
Radius = 500; % Radius of the cells (in m)

p = Radius^(3.8)*10^(SNR/10); % power of the pilots for the desired SNR at the cell edge 

M = 100; % number of antennas at the BS
K = 1; % Number of users per BS
N = 8;
Radius = 500; % Radius of the cells (in m)
nrBS = 2; % Number of BS
beamform = 1; % if beamform = 0, w = [1; 1;], i.e., there is no beamforming at the user
delta = 1;
% generate users (SystemPlot)
% generate one tier (7 BS) with one user per BS. The radius of the BS is
% 500 m
Distances = [300 600; 600 300]; % vector of distances for the UEs

betas = 1./(Distances.^(3.8)); % loss factor
sizeBeta = size(betas);

angularSpread = 10; % 10�
%theta2 = rand(6,1)*pi;
%save('thetaVect.mat','theta2');
theta2 = load('thetaVect.mat');

%% Creating the Rb matrix for the UE1
% Rb matrix from user 1 to BS 1
theta = pi/3; 
Rb(:,:,1) = functionOneRingModel(M,angularSpread,theta);
% Rb matrix from user 1 to BS 2
Rb(:,:,3) = functionOneRingModel(M,angularSpread,theta2.theta2(1));
%Ru matrix for user 1 vs BS 1
Ru(:,:,1) = functionOneRingModel(N,angularSpread,theta2.theta2(2));
% Ru matrix from user 1 vs BS 2
Ru(:,:,3) = functionOneRingModel(N,angularSpread,theta2.theta2(3));

 % Rb matrix from user 2 to BS 2
Rb(:,:,4) = functionOneRingModel(M,angularSpread,theta2.theta2(4));
%Ru matrix of user 2 vs BS 1
Ru(:,:,2) = functionOneRingModel(N,angularSpread,theta2.theta2(5));
% Ru matrix from user 2 vs BS 2
Ru(:,:,4) = functionOneRingModel(N,angularSpread,theta2.theta2(6));
%%

theta21 = linspace(0,pi,250);% Vector of AoA for the second user

RbTot = zeros(M,M,Nrealizations);
RnbTot = zeros(M,M,Nrealizations);
RoptTot = zeros(M,M,Nrealizations);
TMSEnb = zeros(Nrealizations,length(theta21));
TMSEb = zeros(Nrealizations,length(theta21));
TMSEopt = zeros(Nrealizations,length(theta21));
for na = 1:length(theta21)
    
    for r = 1:Nrealizations
        % Rb matrix from user 2 to BS 1
        Rb(:,:,2) = functionOneRingModel(M,angularSpread,theta21(na));
       


        for i=1:2
            Rusum(:,:,i) = zeros(N,N);
            for t = 1:2
                Rusum(:,:,i) = Rusum(:,:,i) + p*betas(i)*Ru(:,:,(t-1)*2 + i); % using betas to model the path loss              
            end

            B(:,:,i) = Rusum(:,:,i) + delta*eye(N);
            [V,D] = eig(B(:,:,i));
            Bsqrt(:,:,i) = V*sqrt(D)*ctranspose(V);
            BusqrtInv(:,:,i) = inv(Bsqrt(:,:,i));
        end

        for n = 1:2
            [eigenVect, eigenVal] = eig(Ru(:,:,(n-1)*2 + n));
            wb(:,n) = eigenVect(:,end);
            wnb(:,n) = ones(N,1)/sqrt(N);

            [V,D] = eig(BusqrtInv(:,:,n)*Ru(:,:,(n-1)*2+n)*BusqrtInv(:,:,n));
            [m,I] = max(abs(diag(D)));
            wopt_ = V(:,I);
            wopt(:,n) = BusqrtInv(:,:,n)*wopt_;
            wopt(:,n) = wopt(:,n)/norm(wopt(:,n));
        end

        for t = 1:2
            for u = 1:2
                Rk_b(:,:,(t-1)*K*nrBS+u) = Ru(:,:,(t-1)*K*nrBS+u)*wb(:,u)*ctranspose(wb(:,u));
                Rk_nb(:,:,(t-1)*K*nrBS+u) = Ru(:,:,(t-1)*K*nrBS+u)*wnb(:,u)*ctranspose(wnb(:,u));
                Rk_opt(:,:,(t-1)*K*nrBS+u) = Ru(:,:,(t-1)*K*nrBS+u)*wopt(:,u)*ctranspose(wopt(:,u));
                %Rk_optCHOL(:,:,(t-1)*K*nrBS+u) = Ru(:,:,(t-1)*K*nrBS+u)*wopt2(:,u)*ctranspose(wopt2(:,u));

                Rkkb(:,:,(t-1)*K*nrBS+u) = betas(u,t)*Rb(:,:,(t-1)*K*nrBS+u)*trace(Rk_b(:,:,(t-1)*K*nrBS+u)); % betas to model the path loss
                Rkknb(:,:,(t-1)*K*nrBS+u) = betas(u,t)*Rb(:,:,(t-1)*K*nrBS+u)*trace(Rk_nb(:,:,(t-1)*K*nrBS+u));
                Rkkopt(:,:,(t-1)*K*nrBS+u) = betas(u,t)*Rb(:,:,(t-1)*K*nrBS+u)*trace(Rk_opt(:,:,(t-1)*K*nrBS+u));
            end
        end

        for t=1:nrBS
            Rsumb = sum(Rkkb(:,:,(t-1)*K*nrBS+1:(t-1)*K*nrBS+nrBS*K),3);
            Rsumnb = sum(Rkknb(:,:,(t-1)*K*nrBS+1:(t-1)*K*nrBS+nrBS*K),3);
            Rsumopt = sum(Rkkopt(:,:,(t-1)*K*nrBS+1:(t-1)*K*nrBS+nrBS*K),3);

            RbTot(:,:,r) = RbTot(:,:,r) + Rkkb(:,:,(t-1)*nrBS + t);
            RnbTot(:,:,r) = RnbTot(:,:,r) + Rkknb(:,:,(t-1)*nrBS + t);
            RoptTot(:,:,r) = RoptTot(:,:,r) + Rkkopt(:,:,(t-1)*nrBS + t);

            for a=1:K
                % For the total NMSE
                Cb(:,:,t,a,r) = Rkkb(:,:,(t-1)*K*nrBS+K*(t-1)+a) - p*Rkkb(:,:,(t-1)*K*nrBS+K*(t-1)+a)/(p*Rsumb + eye(M))*Rkkb(:,:,(t-1)*K*nrBS+K*(t-1)+a);
                Cnb(:,:,t,a,r) = Rkknb(:,:,(t-1)*K*nrBS+K*(t-1)+a) - p*Rkknb(:,:,(t-1)*K*nrBS+K*(t-1)+a)/(p*Rsumnb + eye(M))*Rkknb(:,:,(t-1)*K*nrBS+K*(t-1)+a);
                Copt(:,:,t,a,r) = Rkkopt(:,:,(t-1)*K*nrBS+K*(t-1)+a) - p*Rkkopt(:,:,(t-1)*K*nrBS+K*(t-1)+a)/(p*Rsumopt + eye(M))*Rkkopt(:,:,(t-1)*K*nrBS+K*(t-1)+a);

                TMSEnb(r,na) = TMSEnb(r,na) + trace(Cnb(:,:,t,a,r));
                TMSEb(r,na) = TMSEb(r,na) + trace(Cb(:,:,t,a,r));
                TMSEopt(r,na) = TMSEopt(r,na) + trace(Copt(:,:,t,a,r));

                % For each BS, CALCULATE THE NMSE
                normFactorb = trace(Rkkb(:,:,(t-1)*K*nrBS+K*(t-1)+a));
                normFactornb = trace(Rkknb(:,:,(t-1)*K*nrBS+K*(t-1)+a));
                normFactoropt = trace(Rkkopt(:,:,(t-1)*K*nrBS+K*(t-1)+a));


                MSEb(t,r,a,na) = trace(Cb(:,:,t,a,r))/normFactorb;
                MSEnb(t,r,a,na) = trace(Cnb(:,:,t,a,r))/normFactornb;
                MSEopt(t,r,a,na) = trace(Copt(:,:,t,a,r))/normFactoropt;
            end
        end
        normFactorbTot = trace(RbTot(:,:,r));
        normFactornbTot = trace(RnbTot(:,:,r));
        normFactoroptTot = trace(RoptTot(:,:,r));
        
        TMSEnb(r,na) = TMSEnb(r,na)/normFactornbTot;
        TMSEb(r,na) = TMSEb(r,na)/normFactorbTot;
        TMSEopt(r,na) = TMSEopt(r,na)/normFactoroptTot;
    end
    
    meanMSEb(:,na,a) = mean(MSEb(:,:,1,na),2);
    meanMSEnb(:,na,a) = mean(MSEnb(:,:,1,na),2);
    meanMSEopt(:,na,a) = mean(MSEopt(:,:,1,na),2);
    
end
%%
figure;
hold on
grid on
plot((theta21)*180/pi,10*log10(real(meanMSEnb(1,:,1))),'-o');
plot((theta21)*180/pi,10*log10(real(meanMSEb(1,:,1))),'-+');
plot((theta21)*180/pi,10*log10(real(meanMSEopt(1,:,1))),'-*');


legend('Without Beamforming', 'Max. SNR beamforming', 'Max. SLNR beamforming')

xlabel('AoA of UE 2 IN BS 1')
ylabel('MSE(dB)')
title(['MSE AT BS 1 WITH ' num2str(N) ' ANTENNA AT UE']);
%% 
%%
%%

load('D4.mat');
%%

figure;
hold on
grid on
plot((theta21)*180/pi,10*log10(real(meanMSEnb(1,:,1))),'-o');
plot((theta21)*180/pi,10*log10(real(meanMSEb(1,:,1))),'-+');
plot((theta21)*180/pi,10*log10(real(meanMSEopt(1,:,1))),'-*');

plot((theta21)*180/pi,10*log10(real(mnb4)),'-o');
plot((theta21)*180/pi,10*log10(real(mb4)),'-+');
plot((theta21)*180/pi,10*log10(real(mopt4)),'-*');

legend('Without Beamforming', 'Max. SNR beamforming', 'Max. SLNR beamforming')

xlabel('AoA of UE 2 IN BS 1')
ylabel('MSE(dB)')
title(['MSE AT BS 1 WITH ' num2str(N) ' ANTENNA AT UE']);
%%


mnb4 = meanMSEnb(1,:,1);
mb4 = meanMSEb(1,:,1);
mopt4 = meanMSEopt(1,:,1);
%%
%%
%%
TMSEnb2 = mean(TMSEnb);
TMSEb2 = mean(TMSEb);
TMSEopt2 = mean(TMSEopt);

%%
figure;
grid on
hold on
plot((theta21)*180/pi,10*log10(real(TMSEnb2)));
plot((theta21)*180/pi,10*log10(real(TMSEb2)));
plot((theta21)*180/pi,10*log10(real(TMSEopt2)));

legend('Without Beamforming', 'Max. SNR beamforming', 'Max. SLNR beamforming')

xlabel('AoA difference (degrees)')
ylabel('MSE(dB)')
title(['TOTAL NMSE OF THE NETWORK']);
%%

