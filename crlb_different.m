  clear all
%% �������
M=64;%ofdm���Ÿ���
N=512;%���ز�����
c=3*1e8;%����
delta_f1=30*1e3;%��Ƶ���ز����
delta_f2=120*1e3;%��Ƶ���ز����
fc1=5.9*1e9;%��Ƶ�ز�Ƶ��
fc2=24*1e9;%��Ƶ�ز�Ƶ��
T1=1/delta_f1+ 5.85310734463*10^(-6);%��Ƶ�����ܳ���
T2=1/delta_f2+ 1.3*10^(-6);%��Ƶ�����ܳ���
K=4;%��״��Ƶ���
Q=4;%��״��Ƶ���
Nr=N/K;
MB=M/Q;
snr=-30:-20;
SNR=10.^(snr/10); 
%% CRLB_R1 and CRLB_v1  �߿����
crlb_r1=(3*c^2./(8*pi^2*delta_f2^2.*SNR))*(1/(M*Nr*(Nr-1)*(2*Nr-1)+MB*N*(N-1)*(2*N-1)-(9*(Nr*(Nr-1)*M*(M-1)+N*(N-1)*MB*(MB-1)*Q)^2/(4*(M*Nr*(M-1)*(2*M-1)+N*Q^2*MB*(MB-1)*(2*MB-1))))));
crlb_v1=(3*c^2./(8*pi^2*fc1^2*T1^2.*SNR))*(1/(Nr*M*(M-1)*(2*M-1)+N*Q^2*MB*(MB-1)*(2*MB-1)-(9*(Nr*(Nr-1)*M*(M-1)+N*(N-1)*MB*(MB-1)*Q)^2/(4*(M*Nr*(Nr-1)*(2*Nr-1)+MB*N*(N-1)*(2*N-1))))));
%% CRLB_R2 and CRLB_V2  ����Ϳ�
crlb_r2=(3*c^2./(8*pi^2.*SNR))*(1/(delta_f2^2*K^2*M*Nr*(Nr-1)*(2*Nr-1)+delta_f1^2*MB*N*(N-1)*(2*N-1)-(9*(delta_f2*K*Nr*(Nr-1)*M*(M-1)+delta_f1*N*(N-1)*MB*(MB-1)*Q)^2/(4*(Nr*M*(M-1)*(2*M-1)+N*Q^2*MB*(MB-1)*(2*MB-1))))));
crlb_v2=(3*c^2./(8*pi^2*fc1^2*T1^2.*SNR))*(1/(Nr*M*(M-1)*(2*M-1)+N*Q^2*MB*(MB-1)*(2*MB-1)-(9*(delta_f2*K*Nr*(Nr-1)*M*(M-1)+delta_f1*N*(N-1)*MB*(MB-1)*Q)^2/(4*(delta_f2^2*K^2*M*Nr*(Nr-1)*(2*Nr-1)+delta_f1^2*MB*N*(N-1)*(2*N-1))))));
%% CRLB_R3 and CRLB_V3 �ߵͿ�
crlb_r3=(3*c^2./(8*pi^2.*SNR))*(1/(N*(N-1)*MB*((2*N-1)*(delta_f1^2+delta_f2^2)-(9*(delta_f1+delta_f2)^2*(N-1)*(MB-1)/((8)*(2*MB-1))))));
crlb_v3=(3*c^2./(8*pi^2*fc1^2*T1^2.*SNR))*(1/(N*MB*(MB-1)*(2*(2*MB-1)-(9*(delta_f1+delta_f2)^2*(N-1)*(MB-1)/(4*(delta_f1^2+delta_f2^2)*(2*N-1))))));
%% CRLB_R4 and CRLB_V4  �ߵ���
crlb_r4=(3*c^2./(8*pi^2*delta_f2^2.*SNR))*(1/(Nr*(Nr-1)*M*((2*Nr-1)*(K^2+1)-(9*(1+K)^2*(Nr-1)*(M-1)/((8)*(2*M-1))))));
crlb_v4=(3*c^2./(8*pi^2*fc1^2*T1^2.*SNR))*(1/(Nr*(M-1)*M*(2*(2*M-1)-(9*(1+K)^2*(Nr-1)*(M-1)/(4*(1+K^2)*(2*Nr-1))))));
%% ��Ƶ��״
crlb_lowshu_r=(3*c^2./(8*pi^2*delta_f2^2.*SNR))*(1/(Nr*(Nr-1)*M*(2*Nr-1-((9*(Nr-1)*(M-1))/(4*(2*M-1))))));
crlb_lowshu_v=(3*c^2./(8*pi^2*fc1^2*T1^2.*SNR))*(1/(M*(M-1)*Nr*(2*M-1-((9*(Nr-1)*(M-1))/(4*(2*Nr-1))))));
%% ��Ƶ��״
crlb_highshu_r=(3*c^2./(8*pi^2*delta_f2^2*K^2.*SNR))*(1/(Nr*(Nr-1)*M*(2*Nr-1-((9*(Nr-1)*(M-1))/(4*(2*M-1))))));
crlb_highshu_v=(3*c^2./(8*pi^2*fc1^2*T1^2.*SNR))*(1/(M*(M-1)*Nr*(2*M-1-((9*(Nr-1)*(M-1))/(4*(2*Nr-1))))));
%% ��Ƶ��״
crlb_lowkuai_r=(3*c^2./(8*pi^2*delta_f1^2.*SNR))*(1/(N*(N-1)*MB*(2*N-1-((9*(N-1)*(MB-1))/(4*(2*MB-1))))));
crlb_lowkuai_v=(3*c^2./(8*pi^2*fc1^2*T1^2*Q^2.*SNR))*(1/(MB*(MB-1)*N*(2*MB-1-((9*(N-1)*(MB-1))/(4*(2*N-1))))));
%% ��Ƶ��״
crlb_highkuai_r=(3*c^2./(8*pi^2*delta_f2^2.*SNR))*(1/(N*(N-1)*MB*(2*N-1-((9*(N-1)*(MB-1))/(4*(2*MB-1))))));
crlb_highkuai_v=(3*c^2./(8*pi^2*fc1^2*T1^2*Q^2.*SNR))*(1/(MB*(MB-1)*N*(2*MB-1-((9*(N-1)*(MB-1))/(4*(2*N-1))))));
% %% �Ա�
% figure(1)%��Ƶ��״���Ƶ��״�Ա�
% plot(snr,crlb_v1,'blue')
% hold on
% plot(snr,crlb_v4,'red')
% legend('gaokuai','gaoshu')
%% ����
figure(1)%range
plot(snr,crlb_r1.^0.5,"blue");
hold on
plot(snr,crlb_r2.^0.5,"red");
hold on
plot(snr,crlb_r3.^0.5,"yellow");
hold on
plot(snr,crlb_r4.^0.5,"green");
hold on
xlabel('SNR')
ylabel('CRLB of range')
legend('CA-based staggered pilot','high-comb,low-block','high and low full block','high and low full comb')
figure(2)
plot(snr,crlb_v1.^0.5,"blue");
hold on
plot(snr,crlb_v2.^0.5,"red");
hold on
plot(snr,crlb_v3.^0.5,"yellow");
hold on
plot(snr,crlb_v4.^0.5,"green");
hold on