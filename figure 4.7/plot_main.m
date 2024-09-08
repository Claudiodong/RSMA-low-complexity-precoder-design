clc
clear
close all;
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels

SNR_dB = 5:5:40;
load("2users_4Nr_24Nt_100.mat")

figure()
plot(SNR_dB,Avg_RSMA_1,'--o','LineWidth',1,'MarkerSize',8,'Color',"#D95319");
hold on;
plot(SNR_dB,Avg_RSMA_2,'--x','LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB,Avg_RSMA_3,'--V','LineWidth',1,'MarkerSize',8,'Color',"#EDB120");
hold on;

plot(SNR_dB,Avg_RSMA_4,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#77AC30");


load("3users_4Nr_24Nt_100.mat")


plot(SNR_dB,Avg_RSMA_1,'--o','LineWidth',1,'MarkerSize',8,'Color',"#D95319");
hold on;
plot(SNR_dB,Avg_RSMA_2,'--x','LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB,Avg_RSMA_3,'--V','LineWidth',1,'MarkerSize',8,'Color',"#EDB120");
hold on;
plot(SNR_dB,Avg_RSMA_4,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#77AC30");


load("4users_4Nr_24Nt_100.mat")

plot(SNR_dB,Avg_RSMA_1,'--o','LineWidth',1,'MarkerSize',8,'Color',"#D95319");
hold on;
plot(SNR_dB,Avg_RSMA_2,'--x','LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB,Avg_RSMA_3,'--V','LineWidth',1,'MarkerSize',8,'Color',"#EDB120");
hold on;
plot(SNR_dB,Avg_RSMA_4,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#77AC30");
grid on;
xlabel("SNR [dB]");ylabel("Sum rate [bps/Hz]")

legend("RSMA low complexity $Q_{c}=1$",'RSMA low complexity $Q_{c}=2$','RSMA low complexity $Q_{c}=3$','RSMA low complexity $Q_{c}=4$');