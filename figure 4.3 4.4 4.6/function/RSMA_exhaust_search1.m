function [optimal_RSMA,optimal_t_by_search]=RSMA_exhaust_search1(H,Num_U,P_k,SNR,Rnn,Q_c)
    t = 0:0.001:1;
    R_sum_RSMA = zeros(1,length(t));
    for i_t = 1:length(t)
         [R_sum_RSMA(i_t),~,~]=adaptive_common_stream(H,Q_c,P_k,SNR,t(i_t),Rnn,Num_U);
    end
    
    [optimal_RSMA,index] = max(R_sum_RSMA);
    optimal_t_by_search = t(index);   
end