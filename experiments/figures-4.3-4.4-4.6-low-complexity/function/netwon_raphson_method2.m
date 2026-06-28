function [t] = netwon_raphson_method2(H,P_k,Num_U,tolerance,SNR,Q_c)
    [~,~,V]=svd(H);
    [Nr,~] = size(H);
    N_k = Nr/Num_U;
    P_c = V(:,1:Q_c);
    % netwon raphson method
    for i = 1:Num_U
        H_k = H(N_k*i-(N_k-1):N_k*i,:);

        % private part
        A_term = H_k * P_k(:,:,i);
        A_termss(:,:,i) = A_term* A_term';
    
        % common part
        B_term = H_k * P_c;
        B_terms(:,:,i) = B_term*B_term';
        SINR_k(i) = real( det( B_terms(:,:,i)/(A_termss(:,:,i) + eye(N_k)) ) ) ;
    end
    
    % predetermine the user index k'
    [~,B_index] = min(SINR_k);
    B = B_terms(:,:,B_index); % common 
    A = A_termss(:,:,B_index);% private
    t = 0.001;
    old_t = 0;
    count = 0;
    while(1)
        count = count +1;
        term1 = (SNR/Nr) * A - (SNR/Q_c) * B;
        term2 = eye(N_k) + (SNR*t/Nr) * A + (SNR*(1-t)/Q_c) * B;
%         term2 = (SNR*t/Nr) * A + (SNR*(1-t)/Q_c) * B;

        term5 = 0;
        term6 = 0;
        % rest of private terms
        for i = 1:Num_U
            if (i ~= B_index)
                term3 = (SNR/Nr) * A_termss(:,:,i); 
                term4 = (eye(N_k) + (SNR*t/Nr) * A_termss(:,:,i) );
%                 term4 = (SNR*t/Nr) * A_termss(:,:,i) ;
                term5 = term5 + real(trace(term3 /term4));
                term6 = term6 + real(trace(-(term3*term3) /(term4*term4)));
            end
        end

        first_order = trace(term1/term2) + term5;
%         first_order = trace(term1/term2) + N_k*(Num_U-1)/t;
        second_order = trace(- (term1*term1)/(term2*term2)) + term6;
%         second_order = trace(- (term1*term1)/(term2*term2)) - 2*N_k*(Num_U-1)/(t^(2));
        t = t - first_order/(second_order) ; % adding 1 to ensure the result is at least to be 1
        if(abs(t - old_t)/t < tolerance)
          break
        end
        if(count>1e4)
            break;
        end
        old_t = t;
    end
    % make sure the value of t is between 0 and 1
    t = min(real(t),1);

    % if the t is still less than 0, make it as 0
    if (t <0)
        t = 0;
    end
 
end