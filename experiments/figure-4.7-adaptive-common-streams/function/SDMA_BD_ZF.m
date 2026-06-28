function [R_k_SDMA_BD,Pri_precoder]=SDMA_BD_ZF(SNR,H,P_k,N_k)
        [Nr,Nt]=size(H);
        a_k_SDMA = sqrt(SNR/Nr);
        Num_U = Nr / N_k;

        for i = 1:Num_U
            Pri_precoder(:,N_k*i-(N_k-1):i*N_k) = a_k_SDMA*eye(Nt)*P_k(:,:,i);
        end
 
        % if want to be fair, then should apply same design as the RSMA and
        % SDMA
        R_k_SDMA_BD = zeros(1,Num_U);

        % try to use the zero forcing ?
         for i = 1:Num_U
                B = N_k*i-(N_k-1):N_k*i;
                H_k = H(B,:);
                % BD
                F_BD = a_k_SDMA * eye(N_k)*H_k;
            
                % SNR
                SNR_k_SDMA_BD = (F_BD * P_k(:,:,i)* P_k(:,:,i)'*F_BD');
                R_k_SDMA_BD(i) = real(log2 (det( eye(N_k) + SNR_k_SDMA_BD ))); % Rate
         end
end