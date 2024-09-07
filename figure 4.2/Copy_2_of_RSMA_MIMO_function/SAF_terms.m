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
% - G_MMSE_c_k
% - G_MMSE_p_k (Nr*Nr*Num_U, compelx) = private stream MMSE Combiner design 
% - U_c_k
% - U_p_k (Nr*Nr*Num_U, compelx) = instantaneous weight for WMMSE 
% - sigma (1*NUm_U) = user received noise power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: (SAF terms after vectersiation)
% - A_dot_c_k = (Num_U*Nr*Nt * Num_U*Nr*Nt * Num_U) hermition matrix
% - A_c_k = (Num_U*Nr*Nt * Num_U*Nr*Nt * Num_U) hermition matrix
% - A_p_k (Num_U*Nr*Nt * Num_U*Nr*Nt * Num_U) hermition matrix
% - a_c_k (Num_U*Nr*Nt * Num_U) complex
% - a_p_k (Num_U*Nr*Nt * Num_U) complex
% - phi_c_k (1* Num_U) complex
% - phi_p_k (1* Num_U) complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A_dot_c_k,A_c_k,A_p_k,a_c_k,a_p_k,phi_c_k,phi_p_k]=SAF_terms(H,G_MMSE_c_k,G_MMSE_p_k,U_c_k,U_p_k,sigma,Common_precoder)
          [Nt,~] = size(H);
          [Nr,N_common,Num_U] = size(G_MMSE_c_k);
          N_private = N_common;
          % define the variable
          A_dot_c_k = zeros(Nt*Nr,Nt*Nr,Num_U);
          A_c_k = zeros(Nt*Nr,Nt*Nr,Num_U);
          A_p_k = zeros(Nt*Nr,Nt*Nr,Num_U);
          a_c_k = zeros(Nt*Nr,Num_U);
          a_p_k = zeros(Nt*Nr,Num_U);
          phi_c_k = zeros(1,Num_U);
          phi_p_k = zeros(1,Num_U);

%           power = diag(Common_precoder'*Common_precoder);
%           if(power(1)>0 && power(2)>0)
%               N_common = 2;
%           else
%               N_common = 1;
%           end
    
          % compute the variable
          for i = 1:Num_U
            B = i*2-1:2*i;
            A_dot_c_k1 =kron( eye(N_common) , H(:,B)*G_MMSE_c_k(:,:,i)'*U_c_k(:,:,i)*G_MMSE_c_k(:,:,i)*H(:,B)' ) ;
            A_dot_c_k(:,:,i) = 0.5 * (A_dot_c_k1 + A_dot_c_k1');% make sure it is hermition otherwise unable to process
    
            A_c_k1 = kron( eye(N_private) , H(:,B)*G_MMSE_c_k(:,:,i)'*U_c_k(:,:,i)*G_MMSE_c_k(:,:,i)*H(:,B)' ) ;
            A_c_k(:,:,i) = 0.5 * (A_c_k1 + A_c_k1');% make sure it is hermition
    
            A_p_k1 = kron( eye(N_private) , H(:,B)*G_MMSE_p_k(:,:,i)'*U_p_k(:,:,i)*G_MMSE_p_k(:,:,i)*H(:,B)' ) ;
            A_p_k(:,:,i) = 0.5 * (A_p_k1 + A_p_k1');% make sure it is hermition
    
            a_c_k(:,i) = vec(  H(:,B)  * G_MMSE_c_k(:,:,i)' * U_c_k(:,:,i));
            a_p_k(:,i) = vec(  H(:,B)  * G_MMSE_p_k(:,:,i)' * U_p_k(:,:,i));

    
            phi_c_k(:,i) = real(sigma(i) * trace(U_c_k(:,:,i) * G_MMSE_c_k(:,:,i) * G_MMSE_c_k(:,:,i)' ) + trace(U_c_k(:,:,i)) - log2(det(U_c_k(:,:,i))));
            phi_p_k(:,i) = real(sigma(i) * trace(U_p_k(:,:,i) * G_MMSE_p_k(:,:,i) * G_MMSE_p_k(:,:,i)' ) + trace(U_p_k(:,:,i)) - log2(det(U_p_k(:,:,i))));
          end

end