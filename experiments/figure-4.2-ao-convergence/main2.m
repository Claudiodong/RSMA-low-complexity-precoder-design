clc
clear
close all;
%%
addpath('Copy_2_of_RSMA_MIMO_function\')
addpath('RSMA_MIMO_function\')
% RSMA MIMO BC system model
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%%
% analysis capacity region, under perfect CSIT condition
Nt = 4;                % Transmit antenna
Nr = 2;                % Received antenna per user 
Num_U = 2;             % Number of user
N_common = min(Nt,Nr); % The number of common stream (is 1 when MISO)
N_private = min(Nt,Nr);% number of private stream == Nr
SNR = 10^(5./10);      % SNR level as 20 dB
sigma = [1,1];         % noise power for user 1 & 2
tolerance = 1e-5;      % tolerance value for the iteration
max = 1;               % number of channel realisation
max_count = 1e3;       % maximum number of iteration count

% Define the weight for each user
weight2 = 1;                             % weight for UE2
weight1 = ones(1,length(weight2));       % weight for UE1
weight = [weight1;weight2];              % weight for all
alpha = 0.8;                             % power for private stream

H = (randn(Nt,Num_U*Nr,max) + 1i*randn(Nt,Num_U*Nr,max))./sqrt(2);

%% SNR=5 dB
[R_sum_BD_5]=RSMA_MIMO_rate_BD(H,Nr,SNR,sigma,weight,tolerance,max_count,alpha);
[R_sum_MRT_5]=RSMA_MIMO_rate_convergence(H,Nr,SNR,sigma,weight,tolerance,max_count,alpha);

%% SNR=20 dB
SNR = 10^(20./10);
[R_sum_BD_20]=RSMA_MIMO_rate_BD(H,Nr,SNR,sigma,weight,tolerance,max_count,alpha);
[R_sum_MRT_20]=RSMA_MIMO_rate_convergence(H,Nr,SNR,sigma,weight,tolerance,max_count,alpha);

%% SNR=35 dB
SNR = 10^(35./10);
[R_sum_BD_35]=RSMA_MIMO_rate_BD(H,Nr,SNR,sigma,weight,tolerance,max_count,alpha);
[R_sum_MRT_35]=RSMA_MIMO_rate_convergence(H,Nr,SNR,sigma,weight,tolerance,max_count,alpha);

%% figure plot                          
figure()
plot(1:22,abs(R_sum_MRT_5(1:22)),'r--x')
hold on;
plot(1:22,abs(R_sum_BD_5(1:22)),'b--o')
hold on;
plot(1:22,abs(R_sum_MRT_20(1:22)),'r--x')
hold on;
plot(1:22,abs(R_sum_BD_20(1:22)),'b--o')
hold on;
plot(1:22,abs(R_sum_MRT_35(1:22)),'r--x')
hold on;
plot(1:22,abs(R_sum_BD_35(1:22)),'b--o')
legend('MRT ','BD ','location','northwest');
xlabel('AO iterative')
ylabel('Sum Rate [bps/Hz]')
ylim([0,35])

% save('convexity_data','R_sum_MRT_5','R_sum_BD_5',"R_sum_BD_20","R_sum_MRT_20","R_sum_BD_35","R_sum_MRT_35")

