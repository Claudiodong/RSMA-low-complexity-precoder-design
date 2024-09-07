clc
clear
close all

load("SNR_result_with_SVD_MRT_initilisation.mat")
Avg_RSMA_MRT = Avg_RSMA;
Avg_SDMA_MRT = Avg_SDMA;
clear Avg_SDMA
clear Avg_RSMA
load('SNR_data_MIMO_SDMA_RSMA.mat')
Avg_SDMA_BD = Avg_SDMA;
SNR = 5:5:40;

figure()
plot(SNR,Avg_RSMA_MRT,'--x','LineWidth',1,'MarkerSize',8,'Color',"#0072BD")
hold on;
plot(SNR,Avg_RSMA_max,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#0072BD")
hold on;
plot(SNR,Avg_SDMA_BD,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#D95319")
hold on;
plot(SNR,Avg_SDMA_MRT,'--x','LineWidth',1,'MarkerSize',8,'Color',"#D95319")

legend("RSMA SVD-MRT initialisation",'RSMA SVD-BD/MRT initialisation','SDMA BD/MRT initialisation','SDMA MRT initialisation');
grid on;
% title("Optimised precoder sum rate with different initilisation")
xlabel('SNR [dB]')
ylabel("Sum rate [bps/Hz]")

increased_RSMA = (abs(Avg_RSMA_MRT(8) - Avg_RSMA_max(8))/Avg_RSMA_MRT(8) ) * 100
Dof_RSMA_BD = Avg_RSMA_max(8) / log2(10^(40/10))
Dof_RSMA_MRT = Avg_RSMA_MRT(8) / log2(10^(40/10))
Dof_RSMA_BD - Dof_RSMA_MRT

increased_SDMA = (abs(Avg_SDMA_MRT(8) - Avg_SDMA_BD(8))/Avg_SDMA_MRT(8) ) * 100
Dof_SDMA_BD = Avg_SDMA_BD(8) / log2(10^(40/10))
Dof_SDMA_MRT = Avg_SDMA_MRT(8) / log2(10^(40/10))
Dof_SDMA_BD - Dof_SDMA_MRT