GMMSE = zeros(nrBS,M,length(p),K,realizations);
for m=1:realizations
    if(mod(m,10) == 0)
        m
    end
    for r=1:length(p)% all the realizations (different SNR)

        for t=1:nrBS % for each BS, calculate the MMSE estimator of the channel

            Rsum = sum(R(:,:,(t-1)*K*nrBS+1:(t-1)*K*nrBS+nrBS*K),3);
            
            % ATENTO AL 1 ANTES DE LA r Y DESPUES DE LA m!!
            GMMSE(t,:,r,1,m) = received(t,:,r,m)*R(:,:,(t-1)*K*nrBS+t)*inv(p(r)*Rsum + eye(M));
           
            MSEH(t,r,m) = immse(GMMSE(t,:,r,1,m),h((t-1)*K*nrBS+t,:));
                        
            C(:,:,t) = R(:,:,(t-1)*K*nrBS+t) - p(r)*R(:,:,(t-1)*K*nrBS+t)*inv(p(r)*Rsum + eye(M))*R(:,:,(t-1)*K*nrBS+t);
            
            if(t==5)
               inverse = inv(p(r)*Rsum + eye(M)); 
            end
            
            MSE(t,r) = trace(C(:,:,t));
        
        end

    end

end
%MSEGmean = mean(MSEG,3);
%MSEHmean = mean(MSEH,3);
MSEHmean = mean(MSEH,3);