%% main results
%% With one user per cell (no intercell interference since the users in the
% same cell use orthogonal pilots). We will see the difference between
% the MSE of the system with and without beamforming, and with the optimal
% beamforming. Now we are having in account the path loss.

%% In mainResults5A, the total MSE is calculated for the whole network, not 
% only for each cell. We compare the TNMSE keeping the UE position fixed
% and only changing the correlation matrices at the UE side.

clear all

Nrealizations = 1;
Nmax = 16; % Maximum number of antennas per terminal in the simulation
SNR = -5; % value of the fixed SNR in dB (power of noise = 1)
Radius = 500; % Radius of the cells (in m)

p = Radius^(3.8)*10^(SNR/10); % power of the pilots for the desired SNR at the cell edge 

M = 100; % number of antennas at the BS
K = 1; % Number of users per BS
%N = [1:2:Nmax]; % Number of antennas per User
N = linspace(1,Nmax,Nmax);
Radius = 500; % Radius of the cells (in m)
nrBS = 7; % Number of BS
beamform = 1; % if beamform = 0, w = [1; 1;], i.e., there is no beamforming at the user
delta = 1;
% generate users (SystemPlot)
% generate one tier (7 BS) with one user per BS. The radius of the BS is
% 500 m
Distances = SystemPlot(nrBS,K,Radius);
thetaB = rand(nrBS*K*nrBS,1)*pi;
thetaU = rand(nrBS*K*nrBS,1)*pi;
betas = 1./(Distances.^(3.8)); % loss factor
sizeBeta = size(betas);
%betas = reshape(betas,[sizeBeta(1)*sizeBeta(2),1]);

angularSpread = 10; % 10�
%%

for i=1:nrBS*K*nrBS
    R(:,:,i) = functionOneRingModel(M,angularSpread,thetaB(i));
end

%%
meanMSEb = zeros(nrBS,length(N));
meanMSEnb = zeros(nrBS,length(N));
meanMSEopt = zeros(nrBS,length(N));

meanMSEoptCHOL = zeros(nrBS,length(N));

TMSEnb = zeros(Nrealizations,length(N));
TMSEb = zeros(Nrealizations,length(N));
TMSEopt = zeros(Nrealizations,length(N));
TMSEopt2 = zeros(Nrealizations,length(N));
for na = 1:length(N)% For all the different values of antennas at the user
    wb = zeros(N(na),K*nrBS);
    wnb = zeros(N(na),K*nrBS);
    wopt = zeros(N(na),K*nrBS);
    wopt2 = zeros(N(na),K*nrBS);
    Ru = zeros(N(na),N(na),nrBS*K*nrBS);
    eigenVect = zeros(N(na),N(na),K*nrBS*nrBS);
    eigenVal = zeros(N(na),N(na),K*nrBS*nrBS);
    Rk_b = zeros(N(na),N(na),K*nrBS);
    Rk_nb = zeros(N(na),N(na),K*nrBS);
    Rk_opt = zeros(N(na),N(na),K*nrBS);
    Rk_optCHOL = zeros(N(na),N(na),K*nrBS);
    Rkkb = zeros(M,M,K*nrBS*nrBS);
    Rkknb = zeros(M,M,K*nrBS*nrBS);
    Rkkopt = zeros(M,M,K*nrBS*nrBS);
    RkkoptCHOL = zeros(M,M,K*nrBS*nrBS);
    gEff = zeros(M,K,Nrealizations);
    Rusum = zeros(N(na),N(na),nrBS*K);
    B = zeros(N(na),N(na),nrBS*K);
    B2 = zeros(N(na),N(na),nrBS*K);
    BcholInv = zeros(N(na),N(na),nrBS*K);
    Bsqrt = zeros(N(na),N(na),nrBS*K);
    BusqrtInv = zeros(N(na),N(na),nrBS*K);
    BusqrtInv2 = zeros(N(na),N(na),nrBS*K);
    Bsqrt2 = zeros(N(na),N(na),nrBS*K);

    Rusqrt = zeros(N(na),N(na),nrBS*K*nrBS);
       
    RbTot = zeros(M,M,Nrealizations);
    RnbTot = zeros(M,M,Nrealizations);
    RoptTot = zeros(M,M,Nrealizations);
    
    

    for r = 1:Nrealizations
        for i = 1:nrBS*K*nrBS % Each user has different Ru for each BS
            
            Ru(:,:,i) = functionOneRingModel(N(na),angularSpread,thetaU(i));
        end
      
        
        % For the optimal beamforming
        delta = 1;
        for i=1:nrBS*K
            Rusum(:,:,i) = zeros(N(na),N(na));
            for t = 1:nrBS
                if (t ~= i) 
                    Rusum(:,:,i) = Rusum(:,:,i) + p*betas(i,t)*Ru(:,:,(t-1)*K*nrBS + i); % using betas to model the path loss
                end
                
            end
            B(:,:,i) = Rusum(:,:,i) + delta*eye(N(na));

            [V,D] = eig(B(:,:,i));
            Bsqrt(:,:,i) = V*sqrt(D)*ctranspose(V);
            BusqrtInv(:,:,i) = inv(Bsqrt(:,:,i));
            
        end

    
        for n = 1:nrBS

            for a = 1:K
                [eigenVect(:,:,(n-1)*K+a),eigenVal(:,:,(n-1)*K+a)] = eig(Ru(:,:,(n-1)*K*nrBS+(n-1)*K + a));
                wb(:,(n-1)*K+a) = eigenVect(:,end,(n-1)*K+a);
                wnb(:,(n-1)*K+a) = ones(N(na),1)/sqrt(N(na)); % without beamforming
                
                % For the optimal beamforming
                [V,D] = eig(BusqrtInv(:,:,(n-1)*K+a)*Ru(:,:,(n-1)*K*nrBS+(n-1)*K + a)*BusqrtInv(:,:,(n-1)*K+a));

                [m,I] = max(abs(diag(D)));
                wopt_ = V(:,I);
                wopt(:,(n-1)*K+a) = BusqrtInv(:,:,(n-1)*K+a)*wopt_;
                wopt(:,(n-1)*K+a) = wopt(:,(n-1)*K+a)/norm(wopt(:,(n-1)*K+a));                
                
            end

        end
    
        for t = 1:nrBS

            for u=1:nrBS*K

                Rk_b(:,:,(t-1)*K*nrBS+u) = Ru(:,:,(t-1)*K*nrBS+u)*wb(:,u)*ctranspose(wb(:,u));
                Rk_nb(:,:,(t-1)*K*nrBS+u) = Ru(:,:,(t-1)*K*nrBS+u)*wnb(:,u)*ctranspose(wnb(:,u));
                Rk_opt(:,:,(t-1)*K*nrBS+u) = Ru(:,:,(t-1)*K*nrBS+u)*wopt(:,u)*ctranspose(wopt(:,u));
                
                if(t==1 && u==1)
                   Rucheck(na,r) = trace(Rk_nb(:,:,(t-1)*K*nrBS+u));
                end
                
                Rkkb(:,:,(t-1)*K*nrBS+u) = betas(u,t)*R(:,:,(t-1)*K*nrBS+u)*trace(Rk_b(:,:,(t-1)*K*nrBS+u)); % betas to model the path loss
                Rkknb(:,:,(t-1)*K*nrBS+u) = betas(u,t)*R(:,:,(t-1)*K*nrBS+u)*trace(Rk_nb(:,:,(t-1)*K*nrBS+u));
                Rkkopt(:,:,(t-1)*K*nrBS+u) = betas(u,t)*R(:,:,(t-1)*K*nrBS+u)*trace(Rk_opt(:,:,(t-1)*K*nrBS+u));

            end

        end
       
        % Calculate the MSE of the given realization
        
        
        for t=1:nrBS % for each BS, calculate the MMSE estimator of the channel
            
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
    
    for a = 1:K
        meanMSEb(:,na,a) = mean(MSEb(:,:,a,na),2); % each row is the mean MSE of a BS
                                 % each column is the mean MSE for a
                                 % different number of antennas at the
                                 % users. The third dimension is the user
                                 % in the BS
        meanMSEnb(:,na,a) = mean(MSEnb(:,:,a,na),2);
        meanMSEopt(:,na,a) = mean(MSEopt(:,:,a,na),2);
    end
    
end

TMSEnb = mean(TMSEnb);
TMSEb = mean(TMSEb);
TMSEopt = mean(TMSEopt);
%% 

figure;
grid on
hold on
plot(N,10*log10(real(TMSEnb)));
plot(N,10*log10(real(TMSEb)));
plot(N,10*log10(real(TMSEopt)));

legend('Without Beamforming', 'Max. SNR beamforming', 'Max. SLNR beamforming')

xlabel('N (Antennas at terminals)')
ylabel('MSE(dB)')
title(['TOTAL NMSE OF THE NETWORK']);

%% 
 
%%
Rucheck2 = mean(Rucheck,2);
%% 

%% Plotting the results
for t = 1:nrBS
    for a = 1:K % for each user in the cell
        plotMeanMSE(meanMSEb,meanMSEnb,meanMSEopt,meanMSEoptCHOL,t,a,N,3) 
    end
end
%%
meanMSEbAllBS = mean(meanMSEb,1);
meanMSEnbAllBS = mean(meanMSEnb,1);
meanMSEoptAllBS = mean(meanMSEopt,1);
meanMSEoptAllBSCHOL = mean(meanMSEoptCHOL,1);

figure;
hold on
plot(N,10*log10(real(meanMSEbAllBS)));
plot(N,10*log10(real(meanMSEnbAllBS)));
plot(N,10*log10(real(meanMSEoptAllBS)));
plot(N,10*log10(real(meanMSEoptAllBSCHOL)));
xlabel('N (Antennas at terminals)')
ylabel('MSE(dB)')
title(['global MSE']);
legend('MSE with beamforming', 'MSE without beamforming', 'MSE optimal beamforming')
%%
figure;
hold on
plot(N,10*log10(real(TMSEb(1,:))));
plot(N,10*log10(real(TMSEnb(1,:))));
plot(N,10*log10(real(TMSEopt(1,:))));
plot(N,10*log10(real(TMSEopt2(1,:))));

legend('MSE with beamforming', 'MSE without beamforming', 'MSE optimal beamforming delta = 1/N', 'MSE optimal beamforming delta = 1')
