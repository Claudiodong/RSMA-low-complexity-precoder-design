
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% - Rate_p_k (1*Num_U) = The rate for each user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Rate_p_k,Pri_precoder]=SDMA_MIMO_rate1(H,Nr,max_count,tolerance,sigma,SNR,weight,Pri_precoder)

    % Find the parameter
    [Nt,Num_U_Nr] = size(H);
    N_common = min(Nt,Nr); % The number of common stream (is 1 when MISO)
    N_private = min(Nt,Nr);% number of private stream == Nr
    Num_U = Num_U_Nr /Nr;
    
    old_rate = 0;
    count = 0;
    while(1)
        count = count + 1;
       % get the WMMSE terms
       [G_MMSE_p_k,U_p_k,~]=SDMA_MIMO_terms(H,Nr,sigma,Pri_precoder);

       % compute the Sample average approximation terms (SAF terms)
       [A_p_k,a_p_k,phi_p_k]=SDMA_SAF_terms(H,G_MMSE_p_k,U_p_k,sigma);

       % update the WMMSE , and find the new precoder design with
       % common rate on each user
       [R_sum,Pri_precoder]=SDMA_MIMO_CVX(H,Nr,SNR,weight,A_p_k,a_p_k,phi_p_k);

       % Compute the accuracy
       accuracy = abs((R_sum - old_rate)/R_sum);
       if (accuracy< tolerance)
           break;
       end
       if (count > max_count)
           break
       end
       old_rate = R_sum;
    end
    % compute the prvate rate 
    [~,~,Rate_p_k]=SDMA_MIMO_terms(H,Nr,sigma,Pri_precoder);
end
