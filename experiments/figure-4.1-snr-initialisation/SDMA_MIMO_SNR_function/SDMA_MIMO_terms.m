
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
% - H (Nt*(NUm_U*Nr), complex) = system channel
% - Nr (1*1 real) = number of received antenna on user
% - sigma (1*NUm_U) = user received noise power
% - Pri_precoder ( (NtxNr) * (NtxNr) ) complex = precoder design for the each
%   user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output: (SAF terms after vectersiation)
% - G_MMSE_p_k (Nr*Nr*Num_U, compelx) = private stream MMSE Combiner design 
% - U_p_k (Nr*Nr*Num_U, compelx) = instantaneous weight for WMMSE 
% - Rate_p_k (1*Num_U) real = user private rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [G_MMSE_p_k,U_p_k,Rate_p_k]=SDMA_MIMO_terms(H,Nr,sigma,Pri_precoder)
    % find the parameter
    [Nt,Num_U_Nr] = size(H);
    N_common = min(Nt,Nr); % The number of common stream (is 1 when MISO)
    N_private = min(Nt,Nr);% number of private stream == Nr
    Num_U = Num_U_Nr /Nr;

 % The noise plus interference covariance matrix
    R_p_k = zeros(Nr,Nr,Num_U);
    for j = 1:Num_U

        A = j*2-1:2*j;
        R_p_k(:,:,j) = sigma(j) * eye(Nr);

        for i = 1:Num_U
            B = i*2-1:2*i;
            % private
            if (i~=j)
                R_p_k(:,:,j) = R_p_k(:,:,j) + H(:,A)'*Pri_precoder(:,B)*Pri_precoder(:,B)'*H(:,A);
            end

        end

    end

    % MMSE filter design
    G_MMSE_p_k = zeros(N_private,Nr,Num_U); % should have size of 2x2 in this case for each user

    for i = 1:Num_U 
     % common filter should have size of  N_common*Nr for each user
     % private should have size of N_private*Nr
     B = i*2-1:2*i;
     G_MMSE_p_k(:,:,i) = Pri_precoder(:,B)' * H(:,B) * ( H(:,B)' * ( Pri_precoder(:,B) * Pri_precoder(:,B)') * H(:,B) + R_p_k(:,:,i))^(-1);
    end


    U_p_k  = zeros(Nr,Nr,Num_U);
    Rate_p_k = zeros(1,Nr);

    % MSE matrix
    for i = 1:Num_U
       B = i*2-1:2*i;
      
       MSE_p_k = ( eye(Nr) + Pri_precoder(:,B)' * H(:,B) * (R_p_k(:,:,i))^(-1) * H(:,B)' * Pri_precoder(:,B) )^(-1);
       
       U_p_k(:,:,i) = ( MSE_p_k )^(-1);
      
       Rate_p_k(i) = real( log2(det(U_p_k(:,:,i))));
    end



end