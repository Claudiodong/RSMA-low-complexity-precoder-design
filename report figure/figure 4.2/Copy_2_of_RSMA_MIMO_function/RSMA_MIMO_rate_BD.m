% Claudio Dong 
% Imperial College London]
% 2024/06/06 
% RSMA MSc Project 


% RSMA MIMO BC optimisation Rate function

% Paper Title "Rate-Splitting Multiple Access for Downlink Multiuser MIMO: 
% Precoder Optimization and PHY-Layer Design 
% by Anup Mishra , Yijie Mao , Member, IEEE, Onur Dizdar , Member, IEEE, 
% and Bruno Clerckx , Fellow, IEEE"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% - H (Nt*(NR*Num_U), complex) = System overall channel
% - Nr (1x1 Real) = NUmber of received antenna on user (assumed user have same number of Nr)
% - tolerance (1x1 real) = The accuracy that required to stop the AO
%   iteration
% - sigma (1* Num_U, real) = The received noise power on each user
% - SNR (1x1 real) = The signal to noise ratio in decimal
% - weight (Num_U *1, real) = The weight ratio for each user to compute the
%   weight sum rate
% - max_count (1x1 real) = maximum number of iteration for the AO algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% - R1 (1*1 real) = The user 1 rate
% - R2 (1*1 real) = The user 2 rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_sum]=RSMA_MIMO_rate_BD(H,Nr,SNR,sigma,weight,tolerance,max_count,H_power)

        [Nt,Num_U_Nr] = size(H);
        N_common = min(Nt,Nr); % The number of common stream (is 1 when MISO)
        N_private = min(Nt,Nr);% number of private stream == Nr
        Num_U = Num_U_Nr /Nr;
        a_c = SNR - SNR^(H_power);
        a_k = SNR^(H_power);

    
        [P_k,~]=BD_design(N_private,H);
        for i = 1:Num_U
            P_K1(:,2*i-1:2*i) = P_k(:,:,i);
        end
        % initilised the precoder (MRT-SVD)
        Pri_precoder = sqrt(a_k/Num_U_Nr) .* P_K1; % private precoder design
        % should have size of  Nt*Num_User for each user
    
        [~,~,V] = svd(H);
        u = V(:,1:N_common); % choose the leftest singular vector based on the number of transmit common stream
    
        Common_precoder = sqrt(a_c/N_common).*(u./norm(u)); % common precoder design
 
        

        old_rate = 0;
        count = 0;
        while(1)
            count = count + 1;
           % get the WMMSE terms
           [G_MMSE_c_k,G_MMSE_p_k,U_c_k,U_p_k,~,~]=RSMA_MIMO_terms(H,Nr,sigma,Common_precoder,Pri_precoder);
           % compute the Sample average approximation terms (SAF terms)
           [A_dot_c_k,A_c_k,A_p_k,a_c_k,a_p_k,phi_c_k,phi_p_k]=SAF_terms(H,G_MMSE_c_k,G_MMSE_p_k,U_c_k,U_p_k,sigma,Common_precoder);

           % update the WMMSE , and find the new precoder design with
           % common rate on each user
           [R_sum(count),part_common_rate,Common_precoder,Pri_precoder]=RSMA_MIMO_CVX(H,Nr,SNR,weight,A_dot_c_k,A_c_k,A_p_k,a_c_k,a_p_k,phi_c_k,phi_p_k);

           accuracy = abs((R_sum(count) - old_rate)/R_sum(count));
           if (accuracy< tolerance)
               break;
           end
           if (count > max_count)
               break
           end
           old_rate = R_sum(count);
        end
        % compute the prvate rate 
        [~,~,~,~,Rate_p_k,~]=RSMA_MIMO_terms(H,Nr,sigma,Common_precoder,Pri_precoder);

        R1 = Rate_p_k(1) + (-part_common_rate(1));
        R2 = Rate_p_k(2) + (-part_common_rate(2));
end