function [R_sum]=RSMA_BD_MRC_rate(H,Num_U,G_k,P_k,SNR,t_optimal_solution)

            [Nr,Nt] = size(H);
            N_k = Nr/Num_U;
            Q_c = min(N_k,Nt);

            % power allocation
            a_k = sqrt(SNR*(t_optimal_solution)/Nr);           % private power
            a_c = sqrt(  SNR * (1-t_optimal_solution) /Q_c); % common power

            R_k = zeros(1,Num_U);
            % User private rate calculation
            for i = 1:Num_U
                B = 2*i-1:2*i;
                % BD terms
                F_k = a_k*eye(N_k)*H(B,:);
%                 F_k = a_k*eye(N_k)*G_k(:,:,i)*H(B,:);
            
                % SNR of private , interference is cancelled out due to BD
                % design
                SNR_k = (F_k * P_k(:,:,i)* P_k(:,:,i)'*F_k') * inv(eye(N_k));
                R_k(i) = real(log2 (det( eye(N_k) + SNR_k ))); % rate
                
            end
            
            %% Common Part

            % Common precoder design
            [~,~,V]=svd(H);
            P_c =  V(:,1:Q_c);

            R_k_c = zeros(1,Num_U);
            % Common combiner design
            for i = 1:Num_U
                H_k = H(2*i-1:2*i,:); % user channel
% 
%                 r_k_c = H_k * P_c; 

                % common combiner design (MRC)
%                 G_k_c = (G_k(:,:,i)')^(-1) * (r_k_c ./ (norm(r_k_c)));

                % compute the common part power 
%                 Comm_power_sqrt = a_c*eye(N_k) *G_k_c' * G_k(:,:,i) * r_k_c;
%                 Comm_power_sqrt = a_c*eye(N_k) * G_k(:,:,i) * H_k * P_c;
                Comm_power_sqrt = a_c*eye(N_k) * H_k * P_c;
                Comm_power = Comm_power_sqrt*Comm_power_sqrt';

                % find the interference from the intended private message to
                % the user kth.
%                 interference_sqrt = a_k*eye(N_k) * G_k_c' *G_k(:,:,i) * H_k * P_k(:,:,i);
%                 interference_sqrt = a_k*eye(N_k) * G_k(:,:,i) * H_k * P_k(:,:,i);
                interference_sqrt = a_k*eye(N_k) * H_k * P_k(:,:,i);
                interference_power = interference_sqrt * interference_sqrt';

                % Compute the SINR
                SINR_k_c = Comm_power / (interference_power + eye(N_k));

                % compute the common rate
                R_k_c(i) = real( log2(det(eye(N_k)+ SINR_k_c)) );
            end

            % Find the minimum rate that require to transmit among the users.
            R_c = min(R_k_c);
            R_sum = R_c + sum(R_k); % sum rate
end