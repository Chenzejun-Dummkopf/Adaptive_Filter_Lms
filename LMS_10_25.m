 % Title : LMS algorithm
% Application : system identification
% Reference : 
clc; % ��ȥһҳ����ڣ�������Ļ���Ͻ�
clear all; % �ӹ����ռ�������б����ͺ���
close all; % �ر����д���
% ���������������в�������������������
K = 100; % �������д���
L = 8; % �˲�������
N = 1000; % ��������
n = 1:N; % ϵͳ��������������
% ������������ϵͳ��ʼ������������������
unknown_w = [0.7753 -0.5286 0.6638 -0.3693 0.6873 -0.0742 0.2465 -0.4815]; % δ֪ϵͳȨ����
%unknown_w = [0.0753 -0.0286 0.0638 -0.0693 0.0873 -0.0742 0.0465 -0.0815]; % δ֪ϵͳȨ����
mu = 0.04; % LMS�̶�����
%ָ���ʼ��
MSE = zeros(1,N); % �������
EMSE = zeros(1,N); % ����������
MSD = zeros(1,N); % ����Ȩֵƫ��
%�������������㷨����������������
for k = 1:K
    % �������������㷨��ʼ������������������
    w = zeros(L,1); % Ȩֵ��ʼ��
    % �����������������źš���������������
    u = randn(1,N); % �����ź�
    %�����������������źš���������������
    d = filter(unknown_w,1,u); % �����ź�
    D = awgn(d,30); % ����30db��˹������
%     D = d + sqrt(1e-3)*randn(1,N); % ����30db��˹������
    %���������������¹��̡�����������
    for i = L:N           
        x = u(i:-1:i-L+1)'; % ��������
        e = D(i) - x'*w; % ����ź�
        ee = d(i) - x'*w; % �������
        if rem(i,2)==0
           w = w + mu*2.2*x*e; % ����Ȩֵ�����㷨
        end
%         if rem(i,3)==0
%            w = w + mu/2*x*e; % Ȩֵ����
%         end
%         if rem(i,5)==0
%            w = w + mu/3*x*e; % Ȩֵ����
%         end
        MSE(i) = MSE(i) + e^2;
        EMSE(i) = EMSE(i) + ee^2;
        MSD(i) = MSD(i) + (norm(w' - unknown_w)/norm(unknown_w))^2;
    end
%     for i = N/4:N           
%         x = u(i:-1:i-L+1)'; % ��������
%         e = D(i) - x'*w; % ����ź�
%         ee = d(i) - x'*w; % �������
% %         if rem(i,4)==0
%            w = w + mu/2*x*e; % Ȩֵ����
% %         end
% %         if rem(i,3)==0
% %            w = w + mu/2*x*e; % Ȩֵ����
% %         end
% %         if rem(i,5)==0
% %            w = w + mu/3*x*e; % Ȩֵ����
% %         end
%         MSE(i) = MSE(i) + e^2;
%         EMSE(i) = EMSE(i) + ee^2;
%         MSD(i) = MSD(i) + (norm(w' - unknown_w)/norm(unknown_w))^2;
%     end
end
%����������������������������������
%������������MSE����������������
MSE = MSE/K;
figure(1);
%subplot(221)
plot(n,10*log10(MSE),'r');
xlabel('��������');
ylabel('MSE(dB)');
legend('LMS');
title('MSE(dB)');
grid on;
%������������EMSE����������������
EMSE = EMSE/K;
figure(2);
%subplot(222)
plot(n,10*log10(EMSE));
xlabel('��������');
ylabel('EMSE(dB)');
legend('LMS');
title('EMSE(dB)');
grid on;
%������������MSD����������������
MSD = MSD/K;
figure(3);
%subplot(223)
plot(n,10*log10(MSD));
xlabel('��������');
ylabel('MSD(dB)');
legend('LMS');
title('MSD(dB)');
grid on;
%������������w����������������
figure(4);
%subplot(224)
plot(unknown_w,'k+');
hold on;
plot(w,'r*');
xlabel('�˲�������L');
ylabel('Ȩֵw');
legend('Actual weights','Estimated weights')
title('Comparison of the actual weights and the estimated weights') ;
axis([0,9,-1,1]);