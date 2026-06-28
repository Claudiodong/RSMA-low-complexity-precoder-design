clc
clear
close all;

load("convexity_data.mat")
index = 1:18;
figure()
plot(index,abs(R_sum_MRT_5(index)),'r--x')
hold on;
plot(index,abs(R_sum_BD_5(index)),'b--o')
hold on;
plot(index,abs(R_sum_MRT_20(index)),'r--x')
hold on;
plot(index,abs(R_sum_BD_20(index)),'b--o')
hold on;
plot(index,abs(R_sum_MRT_35(index)),'r--x')
hold on;
plot(index,abs(R_sum_BD_35(index)),'b--o')
legend('SVD-MRT ','SVD-BD ','location','northwest');
xlabel('AO iterative')
ylabel('Sum Rate [bps/Hz]')
ylim([0,50])