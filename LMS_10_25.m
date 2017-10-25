 % Title : LMS algorithm
% Application : system identification
% Reference : 
clc; % 擦去一页命令窗口，光标回屏幕左上角
clear all; % 从工作空间清除所有变量和函数
close all; % 关闭所有窗口
% ………………运行参数……………………
K = 100; % 独立运行次数
L = 8; % 滤波器阶数
N = 1000; % 采样点数
n = 1:N; % 系统变量，迭代次数
% ………………系统初始化……………………
unknown_w = [0.7753 -0.5286 0.6638 -0.3693 0.6873 -0.0742 0.2465 -0.4815]; % 未知系统权向量
%unknown_w = [0.0753 -0.0286 0.0638 -0.0693 0.0873 -0.0742 0.0465 -0.0815]; % 未知系统权向量
mu = 0.04; % LMS固定步长
%指标初始化
MSE = zeros(1,N); % 均方误差
EMSE = zeros(1,N); % 额外均方误差
MSD = zeros(1,N); % 均方权值偏差
%………………算法……………………
for k = 1:K
    % ………………算法初始化……………………
    w = zeros(L,1); % 权值初始化
    % ………………输入信号……………………
    u = randn(1,N); % 输入信号
    %………………期望信号……………………
    d = filter(unknown_w,1,u); % 期望信号
    D = awgn(d,30); % 加入30db高斯白噪声
%     D = d + sqrt(1e-3)*randn(1,N); % 加入30db高斯白噪声
    %………………更新过程………………
    for i = L:N           
        x = u(i:-1:i-L+1)'; % 输入向量
        e = D(i) - x'*w; % 误差信号
        ee = d(i) - x'*w; % 额外误差
        if rem(i,2)==0
           w = w + mu*2.2*x*e; % 部分权值更新算法
        end
%         if rem(i,3)==0
%            w = w + mu/2*x*e; % 权值更新
%         end
%         if rem(i,5)==0
%            w = w + mu/3*x*e; % 权值更新
%         end
        MSE(i) = MSE(i) + e^2;
        EMSE(i) = EMSE(i) + ee^2;
        MSD(i) = MSD(i) + (norm(w' - unknown_w)/norm(unknown_w))^2;
    end
%     for i = N/4:N           
%         x = u(i:-1:i-L+1)'; % 输入向量
%         e = D(i) - x'*w; % 误差信号
%         ee = d(i) - x'*w; % 额外误差
% %         if rem(i,4)==0
%            w = w + mu/2*x*e; % 权值更新
% %         end
% %         if rem(i,3)==0
% %            w = w + mu/2*x*e; % 权值更新
% %         end
% %         if rem(i,5)==0
% %            w = w + mu/3*x*e; % 权值更新
% %         end
%         MSE(i) = MSE(i) + e^2;
%         EMSE(i) = EMSE(i) + ee^2;
%         MSD(i) = MSD(i) + (norm(w' - unknown_w)/norm(unknown_w))^2;
%     end
end
%………………结果输出……………………
%………………MSE……………………
MSE = MSE/K;
figure(1);
%subplot(221)
plot(n,10*log10(MSE),'r');
xlabel('迭代次数');
ylabel('MSE(dB)');
legend('LMS');
title('MSE(dB)');
grid on;
%………………EMSE……………………
EMSE = EMSE/K;
figure(2);
%subplot(222)
plot(n,10*log10(EMSE));
xlabel('迭代次数');
ylabel('EMSE(dB)');
legend('LMS');
title('EMSE(dB)');
grid on;
%………………MSD……………………
MSD = MSD/K;
figure(3);
%subplot(223)
plot(n,10*log10(MSD));
xlabel('迭代次数');
ylabel('MSD(dB)');
legend('LMS');
title('MSD(dB)');
grid on;
%………………w……………………
figure(4);
%subplot(224)
plot(unknown_w,'k+');
hold on;
plot(w,'r*');
xlabel('滤波器阶数L');
ylabel('权值w');
legend('Actual weights','Estimated weights')
title('Comparison of the actual weights and the estimated weights') ;
axis([0,9,-1,1]);