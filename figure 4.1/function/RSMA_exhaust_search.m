function [optimal_RSMA,optimal_t_by_search]=RSMA_exhaust_search(H,Num_U,G_k,P_k,SNR)
    t = 0:0.001:1;
    R_sum_RSMA = zeros(1,length(t));
    for i_t = 1:length(t)
        [R_sum_RSMA(i_t)]=RSMA_BD_MRC_rate(H,Num_U,G_k,P_k,SNR,t(i_t));
    end
    
    [optimal_RSMA,index] = max(R_sum_RSMA);
    optimal_t_by_search = t(index);   
end