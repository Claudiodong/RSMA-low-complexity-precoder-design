function [R_max_RSMA]=selection_max(SNR,Max,R1_all,R2_all,R3_all,R4_all)
    for i_snr = 1:length(SNR)
        for i_max = 1:Max
            R1 = R1_all(i_max,i_snr);
            R2 = R2_all(i_max,i_snr);
            R3 = R3_all(i_max,i_snr);
            R4 = R4_all(i_max,i_snr);
            R_all = [R1,R2,R3,R4];

            R_max_RSMA(i_max,i_snr) = max(R_all);
        end
    end
end
