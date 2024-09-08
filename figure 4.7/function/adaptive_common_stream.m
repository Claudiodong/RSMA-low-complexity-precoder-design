function [Sum_Rate,P_precoder,w_c]=adaptive_common_stream(H,Q_c,P_k,SNR,t,Rnn,Num_U)

      [Nr,Nt]=size(H);
      % Power design
      a_c = sqrt(SNR*(1-t)/Q_c); % common stream power
      a_k = sqrt(SNR*t/Nr); % private stream power
      N_k = Nr/Num_U;

      % common precoder svd
      [~,~,V] = svd(H); % by using svd 
      P_c = V(:,1:Q_c); % normalised common precoder design based on the number of common stream design
      w_c = a_c .* eye(Nt) * P_c/norm(P_c); % precoder design with power
      R_c_k = zeros(1,Num_U);
      R_k = zeros(1,Num_U);

      % loop for rate calculation
      for i = 1:Num_U
          H_k = H(N_k*i-(N_k-1):N_k*i,:);
          com_sqrt =  H_k * w_c;
          c_power = com_sqrt*com_sqrt'; 

          % compute the interference power term 
          w_k(:,:,i) = a_k*eye(Nt)*P_k(:,:,i);
          priv = H_k * w_k(:,:,i);
          interference = priv*priv';

          %% compute the SNR for single common stream without selection
          SNR_c_k =  c_power / ( interference + Rnn);
          % compute the common rate
          R_c_k(i) = log2( det(Rnn + SNR_c_k ) );

          %% compute the user private rate and SNR
          SNR_k = interference* inv(Rnn);
          R_k(i) = real ( log2( det( Rnn + SNR_k))) ;

          % private precoder for SDMA optimisation
          P_precoder(:,N_k*i-(N_k-1):N_k*i) =  w_k(:,:,i);
      end
      % check the power constraint
      power_constraint_check(w_c,w_k,SNR,Num_U)
      % Sum Rate
      Sum_Rate =real( min(R_c_k) + sum(R_k) );
end