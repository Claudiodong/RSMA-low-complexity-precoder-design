function [t] = netwon_raphson_method(H,P_k,Num_U,tolerance,Q_c)
    [~,~,V]=svd(H);
    P_c = V(:,1:Q_c);
    [Nr,~] = size(H);
    N_k = Nr/Num_U;
    % netwon raphson method
    for i = 1:Num_U
        H_k = H(N_k*i-(N_k-1):N_k*i,:);
        % private part
        A_term = H_k * P_k(:,:,i);
        A_termss(:,:,i) = A_term* A_term';
    
        % commn part
        B_term = H_k * P_c;
        B_terms(:,:,i) = B_term*B_term';
        SINR_k(i) = real( det( B_terms(:,:,i)/(A_termss(:,:,i) + eye(N_k)) ) ) ;
    end
    
    [~,B_index] = min(SINR_k);
    B = B_terms(:,:,B_index); % common 
    A = A_termss(:,:,B_index);% private
    t = 0.0001;
    old_t = 0;
    count = 0;
    while(1)
        count = count +1;
        term1 =  ( (Q_c/Nr)*A - B);
        term2 = ( t * ( (Q_c/Nr) * A - B ) + B ) ;
        first_order = trace(term1 /term2) + Num_U/t;
        second_order = trace( -(term1*term1) / (term2*term2)) - Num_U/((t)^(2));
        t = t - first_order/second_order ;
        if(abs(t - old_t)/t < tolerance)
          break
        end
        if(count>1e5)
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