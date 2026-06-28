function [R_k_SDMA_BD,R_k_SDMA_ZF,Pri_precoder]=SDMA_BD_ZF(SNR,H,P_k,N_k)
        [Nr,Nt]=size(H);
        a_k_SDMA = sqrt(SNR/Nr);
        Num_U = Nr / N_k;

        for i = 1:Num_U
            Pri_precoder(:,2*i-1:2*i) = a_k_SDMA*eye(Nr)*P_k(:,:,i);
        end
 
        % if want to be fair, then should apply same design as the RSMA and
        % SDMA
        R_k_SDMA_ZF = zeros(1,Num_U);
        R_k_SDMA_BD = zeros(1,Num_U);
        P_ZF = zeros(Nr,Nt);

        F_ZF = pinv(H);
        for i = 1:Nt
            P_ZF(:,i) = F_ZF(:,i)/norm(F_ZF(:,i));
        end
        % try to use the zero forcing ?
         for i = 1:Num_U
                B = 2*i-1:2*i;
                H_k = H(B,:);
                % BD
                F_BD = a_k_SDMA * eye(Num_U)*H_k;
            
                % SNR
                SNR_k_SDMA_BD = (F_BD * P_k(:,:,i)* P_k(:,:,i)'*F_BD') * inv(eye(N_k));
                SNR_k_SDMA_ZF = ((a_k_SDMA^2)*eye(N_k)*H_k * P_ZF(:,B) * P_ZF(:,B)' * H_k') * inv(eye(N_k));
                R_k_SDMA_ZF(i) = real(log2 (det( eye(N_k) + SNR_k_SDMA_ZF ))); % Rate
                R_k_SDMA_BD(i) = real(log2 (det( eye(N_k) + SNR_k_SDMA_BD ))); % Rate
         end
end