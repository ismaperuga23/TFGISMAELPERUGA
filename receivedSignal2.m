function y = receivedSignal(p,nrBS,K,M,noisePower,g,N,betas)
    % p: power of the pilot
    % nrBS: number of BS
    % K: Users per BS
    % M: number of antennas at the BS
    % noisePower: power of the noise
    % g: matrix of channel coefficients
    % antennasPerUser: antennas per user
    % betas: loss coefficient for each user
    
    y = zeros(nrBS,M); % The received vector. One row for each BS
    sigma = sqrt(noisePower); %noise standard deviation

    for m=1:nrBS % for each BS, one received signal
  
        for u=1:nrBS*K % total number of user
            noise = (sqrt(2)./2)*(randn(1,M)+1i*randn(1,M)); %noise signal
            sum(g(:,:,(m-1)*K+u),2);
            y(m,:) = y(m,:) + ...
                            sqrt(p)*sqrt(betas(u,m))*(sum(g(:,:,(m-1)*K*nrBS+u),2))' + ...
                            noise; 
%             y(m,:) = y(m,:) + ...
%                             sqrt(p)*sqrt(betas(u,m))*sum(g(:,:,(m-1)*antennasPerUser*nrBS*K+1+(u-1)*antennasPerUser:(m-1)*antennasPerUser*nrBS*K+1+(u-1)*antennasPerUser+(antennasPerUser-1)),2)+ ...
%                             noise; 
        end
        %y(m,:) = sqrt(p)*sqrt(betas(m,m))*sum(g((m-1)*antennasPerUser*nrBS*K+1+(m-1)*antennasPerUser:(m-1)*antennasPerUser*nrBS*K+1+(m-1)*antennasPerUser+(antennasPerUser-1),:),1)+noise;
    end
end