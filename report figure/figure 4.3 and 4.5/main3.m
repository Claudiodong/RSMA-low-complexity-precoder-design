clc
clear
close all
%% Extra path and function
addpath('function\')
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% System Parameter
tic
Nt = 4;                             % Tx number
N_k = 2;                            % Rx number per user
Num_U = 2;                          % number of user
Nr = N_k*Num_U;                     % total number of Rx in the system    
SNR_dB = 5:5:50;                    % SNR in dB
SNR = 10.^(SNR_dB./10);             % snr in linear
Max = 100;                          % number of realisation
tolerance =1e-5;                    % tolerance for iteative algorithm
Rnn = eye(N_k);                     % noise covariance matrix
max_count = 1e3;                    % maximum number of loop for iterative algorithm
%% predefine the variable size
weight = [1;1];
sigma = [1,1];
H_all = (randn(N_k*Num_U,Nt,Max) + 1i* randn(N_k*Num_U,Nt,Max))./sqrt(2);

All_t = zeros(Max,length(SNR_dB));
R_max_SDMA = zeros(Max,length(SNR_dB));
R_1_stream = zeros(Max,length(SNR_dB));
R_2_stream = zeros(Max,length(SNR_dB));
R_stream_max = zeros(Max,length(SNR_dB));
R_SDMA_optimal = zeros(Max,length(SNR_dB));
R_RSMA_optimal = zeros(Max,length(SNR_dB));

% load("H_infor.mat");

%%  Start the loop

for i_snr = 1:length(SNR)
     for i_max = 1:Max
        H = H_all(:,:,i_max);
        % BD precoder design
        [P_k,~]=BD_design(N_k,H); % do not use the combine
      
        %% optimal t by newton raphson method for private part
        [t] = netwon_raphson_method(H,P_k,Num_U,tolerance);% irrevelant to SNR

        %% using newton raphson to determine the number of common stream (min 1, max 2)

        All_t(i_max,i_snr) = t;
        % SDMA Rate
        [R_k_SDMA_BD,~,Pri_precoder_SDMA]=SDMA_BD_ZF(SNR(i_snr),H,P_k,N_k);

        % RSMA Sum rate with different stream
        R_max_SDMA(i_max,i_snr) = sum(R_k_SDMA_BD);

        %  one stream
        [R_1_stream(i_max,i_snr),~,~]=adaptive_common_stream(H,1,Nr,P_k,SNR(i_snr),t,Rnn,Num_U);

        % two streams
        [R_2_stream(i_max,i_snr),Pri_precoder,Com_precoder]=adaptive_common_stream(H,2,Nr,P_k,SNR(i_snr),t,Rnn,Num_U);

        [Optimal_R_2_stream(i_max,i_snr),exhaust_t2(i_max,i_snr)]=RSMA_exhaust_search1(H,Num_U,P_k,SNR(i_snr),Nr,Rnn,2);
        [Optimal_R_1_stream(i_max,i_snr),exhaust_t1(i_max,i_snr)]=RSMA_exhaust_search1(H,Num_U,P_k,SNR(i_snr),Nr,Rnn,1);

        % selection among all the streams, and selec the maximum stream
        % performance
        if (R_1_stream(i_max,i_snr)>R_2_stream(i_max,i_snr))
            R_max_selection(i_max,i_snr) = R_1_stream(i_max,i_snr);
            Q_number(i_max,i_snr) = 1;
        else
            R_max_selection(i_max,i_snr) = R_2_stream(i_max,i_snr);
            Q_number(i_max,i_snr) = 2;
        end
        error(i_max,i_snr) = R_2_stream(i_max,i_snr) - Optimal_R_2_stream(i_max,i_snr);

    end
end 

%% average rate
Avg_R1_RSMA = mean(R_1_stream,1); % proposed with 1 common stream RSMA
Avg_R2_RSMA = mean(R_2_stream,1); % proposed with 2 common stream RSMA
Avg_SDMA_BD = mean(R_max_SDMA,1);
Avg_selective_RSMA = mean(R_max_selection,1);
Avg_RSMA2_exhaust = mean(Optimal_R_2_stream,1); % exhaust search with 2 common stream RSMA
Avg_RSMA1_exhaust = mean(Optimal_R_1_stream,1); % exhaust search with 1 common stream RSMA
Avg_error = mean(mean(error,1));


% %% Compare all result figure 2
% figure()
% plot(SNR_dB,Avg_R1_RSMA,"--x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
% hold on;
% plot(SNR_dB,Avg_R2_RSMA,"--o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
% hold on;
% plot(SNR_dB,Avg_selective_RSMA,"-*",'LineWidth',1,'MarkerSize',8,'Color',"#EDB120"); % select the highest rate among all streams
% hold on;
% plot(SNR_dB,Avg_SDMA_BD,"--v",'LineWidth',1,'MarkerSize',8,'Color',"#7E2F8E"); % SDMA BD rate
% hold on;
% plot(SNR_dB,Avg_RSMA2_exhaust,"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#77AC30");
% legend("RSMA BD with $Q_{c} = 1$","RSMA BD with $Q_{c} = 2$",'RSMA BD with Adaptive $Q_{c}$','SDMA BD','RSMA BD $Q_{c}=2$ with exhaust','Location','northwest');
% grid on;
% xlabel("SNR [dB]")
% ylabel("Sum rate [bps]");

%%
% figure()
% plot(SNR_dB,Avg_R1_RSMA,"--x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
% hold on;
% plot(SNR_dB,Avg_R2_RSMA,"--o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
% hold on;
% plot(SNR_dB,Avg_SDMA_BD,"--v",'LineWidth',1,'MarkerSize',8,'Color',"#7E2F8E"); % SDMA BD rate
% hold on;
% plot(SNR_dB,Avg_RSMA2_exhaust,"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#77AC30");
% legend("RSMA BD with $Q_{c} = 1$","RSMA BD with $Q_{c} = 2$",'SDMA BD','RSMA BD $Q_{c}=2$ with exhaust','Location','northwest');
% grid on;
% xlabel("SNR [dB]")
% ylabel("Sum rate [bps]");

%% figure 4.3 in the report
figure()
plot(SNR_dB,Avg_R1_RSMA,"-x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
hold on;
plot(SNR_dB,Avg_R2_RSMA,"-o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
hold on;
plot(SNR_dB,Avg_RSMA1_exhaust,"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB,Avg_RSMA2_exhaust,"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#D95319");
hold on;
plot(SNR_dB,Avg_SDMA_BD,"--v",'LineWidth',1,'MarkerSize',8,'Color',"#7E2F8E"); % SDMA BD rate
legend("RSMA SVD-BD $Q_{c}=1$ with proposed","RSMA SVD-BD $Q_{c}=2$ with proposed","RSMA SVD-BD $Q_{c}=1$ with exhaust","RSMA SVD-BD $Q_{c}=2$ with exhaust",'SDMA BD','Location','northwest');
grid on;
xlabel("SNR [dB]")
ylabel("Sum rate [bps]");

axes("Position",[0.65,0.15,0.2,0.2]) % location of the small figure and size
box on; % show the box
% area where to magnify
index_x = 5:8;
plot(SNR_dB(index_x),Avg_R1_RSMA(index_x),"-x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
hold on;
plot(SNR_dB(index_x),Avg_R2_RSMA(index_x),"-o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
hold on;
plot(SNR_dB(index_x),Avg_RSMA1_exhaust(index_x),"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB(index_x),Avg_RSMA2_exhaust(index_x),"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#D95319");
axis tight;
grid on;

%% squared error
SE_2 = (Avg_RSMA2_exhaust - Avg_R2_RSMA).^(2);
SE_1 = (Avg_RSMA1_exhaust - Avg_R1_RSMA).^(2);

MSE2 = sum(SE_2)/8;
MSE1 = sum(SE_1)/8;

figure() % squared error
plot(SNR_dB,SE_2,'-x','LineWidth',1,'MarkerSize',8)
hold on;
plot(SNR_dB,SE_1,'-o','LineWidth',1,'MarkerSize',8)
legend("$Q_{c}=2$","$Q_{c}=1$");
grid on;
ylabel("Squared Error");
xlabel("SNR [dB]")

%% compare the selective and proposed adaptive
% figure 4.4 in the report
figure()
max_index = 1:100;
plot(max_index,R_1_stream(max_index,3),'--x','LineWidth',1,'MarkerSize',8);
hold on;
plot(max_index,R_2_stream(max_index,3),'--*','LineWidth',1,'MarkerSize',8);
hold on;
plot(max_index,R_max_selection(max_index,3),'-o','LineWidth',1,'MarkerSize',8);
xlabel("Channel Index");
ylabel("Sum Rate [bps/Hz]");
grid on;
legend("RSMA $Q_{c} = 1$",'RSMA $Q_{c}=2$','RSMA Adaptive $Q_{c}$')

%%
abs(Avg_SDMA_BD - Avg_R2_RSMA)./Avg_SDMA_BD

toc