 %GMMSE(t,:,r) = received(t,:,r)*R(:,:,(t-1)*K*nrBS+t)*inv(p*R(:,:,(t-1)*K*nrBS+t) + eye(M));% without takin in account interferences
            %GMMSE(t,:,r) = received(t,:,r)*R(:,:,(t-1)*K*nrBS+t)*inv(R(:,:,(t-1)*K*nrBS+t)+eye(M));

            %GMMSE(t,:,r) =  received(t,:,r)*(inv(R(:,:,(t-1)*nrBS*antennasPerUser*K+1+(t-1)*antennasPerUser+(a-1))+(noisePower(r)*M*eye(M)))*R(:,:,(t-1)*nrBS*antennasPerUser*K+1+(t-1)*antennasPerUser+(a-1)));
            % calculate the MSE
            % size(MSEH) = nrBS X realizations X antennasPerUser in its
            % cell
            %MSEG(t,r,m) = abs(mean((GMMSE(t,:,r,1,m)-g((t-1)*K*nrBS+t,:)).^2));
            %MSEH(t,r,m) = abs(mean((GMMSE(t,:,r,1,m)-h((t-1)*K*nrBS+t,:)).^2));
          
            
 %         MSEG(r) = abs(mean(mean((GMMSE(:,:,r)-g).^2)));
        %         MSEH(r) = abs(mean(mean((GMMSE(:,:,r)-h).^2)));
