function [P_k]=BD_design(N_k,H)
        [Nr,Nt] = size(H);
        Num_U = Nr/N_k;

        P_k = zeros(Nt,N_k,Num_U);
        G_k = zeros(N_k,N_k,Num_U);
        % do the BD design
        for i = 1:Num_U
            B = 2*i-1:2*i;
            H_k_null = H;
            H_k_null(B,:) = []; % find the interference channel
        
            % do the svd
            [~,~,V_k] = svd(H_k_null);
            V_k_null = V_k(:,N_k+1:end);
        
            H_eff = H(B,:)*V_k_null;
        
            [~,~,V_eff] = svd(H_eff);
        
            V_eff_signal = V_eff(:,1:N_k);
        
            % precoder design
            P_k(:,:,i) = V_k_null*V_eff_signal;
        end
end
