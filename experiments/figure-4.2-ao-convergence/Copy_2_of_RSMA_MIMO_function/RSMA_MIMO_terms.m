
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
% - Common_precoder (Nt* N_COMMON) = common stream precoder design
% - Pri_precoder ( (NtxNr) * (NtxNr) ) complex = precoder design for the each
%   user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output: (SAF terms after vectersiation)
% - G_MMSE_c_k (N_common*N_common*Num_U, compelx) = common stream MMSE Combiner design 
% - G_MMSE_p_k (N_private*N_private*Num_U, compelx) = common stream MMSE Combiner design 
% - U_c_k (Nr*Nr*Num_U, compelx) = instantaneous weight for WMMSE 
% - U_p_k (Nr*Nr*Num_U, compelx) = instantaneous weight for WMMSE 
% - Rate_p_k (1*Num_U) real = user private rate
% - Rate_c (1*Num_U) real = common part rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [G_MMSE_c_k,G_MMSE_p_k,U_c_k,U_p_k,Rate_p_k,Rate_c]=RSMA_MIMO_terms(H,Nr,sigma,Common_precoder,Pri_precoder)


[Nt,Num_U_Nr] = size(H);
N_common = min(Nt,Nr); % The number of common stream (is 1 when MISO)
N_private = min(Nt,Nr);% number of private stream == Nr
Num_U = Num_U_Nr /Nr;

 % The noise plus interference covariance matrix
    R_c_k = zeros(Nr,Nr,Num_U);
    R_p_k = zeros(Nr,Nr,Num_U);
    for j = 1:Num_U

        A = j*2-1:2*j;
        R_c_k(:,:,j) = sigma(j) *eye(Nr); % noise matrix
        R_p_k(:,:,j) = sigma(j) * eye(Nr);

        for i = 1:Num_U
            B = i*2-1:2*i;
            % common 
            R_c_k(:,:,j) = R_c_k(:,:,j) + H(:,A)'*Pri_precoder(:,B)*Pri_precoder(:,B)'*H(:,A);
            % private
            if (i~=j)
                R_p_k(:,:,j) = R_p_k(:,:,j) + H(:,A)'*Pri_precoder(:,B)*Pri_precoder(:,B)'*H(:,A);
            end

        end

    end

    % MMSE filter design
    G_MMSE_c_k = zeros(N_common,Nr,Num_U); % should have size of 2x2 for each user
    G_MMSE_p_k = zeros(N_private,Nr,Num_U); % should have size of 2x2 in this case for each user

    for i = 1:Num_U 
     % common filter should have size of  N_common*Nr for each user
     % private should have size of N_private*Nr
     B = i*2-1:2*i;
     G_MMSE_c_k(:,:,i) = Common_precoder' * H(:,B) * ( H(:,B)' * (Common_precoder*Common_precoder') * H(:,B) + R_c_k(:,:,i) )^(-1);
     G_MMSE_p_k(:,:,i) = Pri_precoder(:,B)' * H(:,B) * ( H(:,B)' * ( Pri_precoder(:,B) * Pri_precoder(:,B)') * H(:,B) + R_p_k(:,:,i))^(-1);
    end

    U_c_k  = zeros(Nr,Nr,Num_U);
    U_p_k  = zeros(Nr,Nr,Num_U);
    Rate_c_k = zeros(1,Nr);
    Rate_p_k = zeros(1,Nr);

    % MSE matrix
    for i = 1:Num_U
       B = i*2-1:2*i;
       MSE_c_k = ( eye(Nr) + Common_precoder'   * H(:,B) * (R_c_k(:,:,i))^(-1) * H(:,B)' * Common_precoder     )^(-1);
       MSE_p_k = ( eye(Nr) + Pri_precoder(:,B)' * H(:,B) * (R_p_k(:,:,i))^(-1) * H(:,B)' * Pri_precoder(:,B) )^(-1);
       U_c_k(:,:,i) = ( MSE_c_k )^(-1);
       U_p_k(:,:,i) = ( MSE_p_k )^(-1);
       Rate_c_k(i) = real( log2(det(U_c_k(:,:,i)))); % choose the minimum common rate such that all the users able to decode it 
       Rate_p_k(i) = real( log2(det(U_p_k(:,:,i))));
    end

    % common rate
    Rate_c = min(Rate_c_k);

end