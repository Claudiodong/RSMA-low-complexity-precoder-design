
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
% - SNR (1x1 real) = The signal to noise ratio in decimal
% - weight (Num_U *1, real) = The weight ratio for each user to compute the
%   weight sum rate
% - A_p_k (Num_U*Nr*Nt * Num_U*Nr*Nt * Num_U) hermition matrix
% - a_p_k (Num_U*Nr*Nt * Num_U) complex
% - phi_p_k (1* Num_U) complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output: (SAF terms after vectersiation)
% - R_sum (1*1 real) = weighted system sum rate of all users
% - p_precoder ((Nt*Nr)*(Nt*Nr)) = private stream precoder design (updated by solve the CVX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_sum,p_precoder]=SDMA_MIMO_CVX(H,Nr,SNR,weight,A_p_k,a_p_k,phi_p_k)

        [Nt,Num_U_Nr] = size(H);
        N_common = min(Nt,Nr); % The number of common stream (is 1 when MISO)
        N_private = min(Nt,Nr);% number of private stream == Nr
        Num_U = Num_U_Nr /Nr;
        % AO optimisation
        cvx_begin quiet
        cvx_solver sedumi
    
          variable p_precoder(Nt,Nr*Num_U) complex;

          % vectorization
          p_k = cvx(zeros(Nt*Nr,Num_U));
          % vectorised the precoder
          for i = 1:Num_U
            B = i*2-1:2*i;
            p_k(:,i) = vec(p_precoder(:,B));
          end
    
          epison_p_k = cvx(zeros(1,Num_U));
            for i = 1:Num_U
                user_interference = 0;
                % compute the user interference
                  for j = 1:Num_U  
                       if (i~=j)
                           user_interference = user_interference + p_k(:,j)'*A_p_k(:,:,i)*p_k(:,j);
                       end
                  end
                  % private epslion
                  part = - a_p_k(:,i)'*p_k(:,i) - p_k(:,i)'*a_p_k(:,i) + phi_p_k(:,i); 
                  epison_p_k(i) = p_k(:,i)'*A_p_k(:,:,i)*p_k(:,i) + user_interference + part;
            end   

            % weighted sum rate
            R_sum =  sum(weight.'.*epison_p_k);
   
            minimise R_sum
            subject to 
              % power constraint
              p_precoder(:)'*p_precoder(:) <= SNR
   
      cvx_end
end


