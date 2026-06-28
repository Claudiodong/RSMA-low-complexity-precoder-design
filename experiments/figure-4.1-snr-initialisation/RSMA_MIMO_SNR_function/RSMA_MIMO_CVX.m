
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
% - SNR (1x1 real) = The signal to noise ratio in decimal
% - weight (Num_U *1, real) = The weight ratio for each user to compute the
%   weight sum rate
% - A_dot_c_k = (Num_U*Nr*Nt * Num_U*Nr*Nt * Num_U) hermition matrix
% - A_c_k = (Num_U*Nr*Nt * Num_U*Nr*Nt * Num_U) hermition matrix
% - A_p_k (Num_U*Nr*Nt * Num_U*Nr*Nt * Num_U) hermition matrix
% - a_c_k (Num_U*Nr*Nt * Num_U) complex
% - a_p_k (Num_U*Nr*Nt * Num_U) complex
% - phi_c_k (1* Num_U) complex
% - phi_p_k (1* Num_U) complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output: (SAF terms after vectersiation)
% - R_sum (1*1 real) = weighted system sum rate of all users
% - part_common_rate (1*Num_U) real = The common rate distributed on each user
% - c_precoder (Nt* N_common) = The common stream precoder design
% - p_precoder ((Nt*N_private)*(Nt*N_private)) = private stream precoder design (updated by solve the CVX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_sum,part_common_rate,c_precoder,p_precoder]=RSMA_MIMO_CVX(H,Nr,SNR,weight,A_dot_c_k,A_c_k,A_p_k,a_c_k,a_p_k,phi_c_k,phi_p_k)
        % AO optimisation
        [Nt,Num_U_Nr] = size(H);
        N_common = min(Nt,Nr); % The number of common stream (is 1 when MISO)
        N_private = min(Nt,Nr);% number of private stream == Nr
        Num_U = Num_U_Nr /Nr;
        cvx_begin quiet
        cvx_solver sedumi
    
          variable c_precoder(Nt,Nr) complex;
          variable p_precoder(Nt,Nr*Num_U) complex;
          variable part_common_rate(1,Num_U);
    
          % vectorization
          p_c = vec(c_precoder);
          p_k = cvx(zeros(Nt*Nr,Num_U));
          for i = 1:Num_U
            B = i*2-1:2*i;
            p_k(:,i) = vec(p_precoder(:,B));
          end
    
          epison_c_k = cvx(zeros(1,Num_U));
          epison_p_k = cvx(zeros(1,Num_U));
            for i = 1:Num_U
                user_interference = 0;
                user_sum = 0;
                  for j = 1:Num_U
                       user_sum = user_sum + p_k(:,j)'*A_c_k(:,:,i) * p_k(:,j);
                       
                       if (i~=j)
                           user_interference = user_interference + p_k(:,j)'*A_p_k(:,:,i)*p_k(:,j);
                       end
                  end
                  A = - a_c_k(:,i)'*p_c - p_c'*a_c_k(:,i)+ phi_c_k(:,i); 
    
                  % common epslion
                  epison_c_k(i) = p_c'*A_dot_c_k(:,:,i)*p_c + user_sum + A;   
    
                  % private epslion
                  part = - a_p_k(:,i)'*p_k(:,i) - p_k(:,i)'*a_p_k(:,i) + phi_p_k(:,i); 
                  epison_p_k(i) = p_k(:,i)'*A_p_k(:,:,i)*p_k(:,i) + user_interference + part;
            end
       
            % weighted sum rate
            R_sum =  sum(weight.'.*(part_common_rate + epison_p_k));
    
            minimise R_sum
            subject to 
              % power constraint
              c_precoder(:)'*c_precoder(:) + p_precoder(:)'*p_precoder(:) <= SNR
              % sum of all common rate at different user should be less
              % than 0 since it is -x.
              for i = 1:Num_U
                sum(part_common_rate) + 2 >= epison_c_k(i)
              end
              part_common_rate <= 0
              
      cvx_end
end


