%% ����˵����������Ƶ�μ��ز��ۺϵĽ���Ƶ����
%���빦�ܣ�
%��д�ˣ�������
%���ڣ�2022.11.30
%% ���������
clear;
clc;
%% ��Ƶ����ϵͳ����
fc1=24e9;	%��Ƶ�ز�Ƶ��
symbols=64; %OFDM����64
carriers=512; %���ز�����
delta_f1=120e3;% ��Ƶ���ز����
c0=3e8;  %����
cp = 1.33*10^(-6);%���������200��
T1=1/delta_f1+cp;%��Ƶ���ų���
%% ��Ƶ����ϵͳ����
fc2=5.9e9;	%��Ƶ�ز�Ƶ��
delta_f2=30e3; %��Ƶ���ز����
T2=1/delta_f2+5.975141242935*10^(-6);%��Ƶ���ų���
%% Ŀ�����
Range_target=117;    %Ŀ�����
Velocity_target=30;  %Ŀ���ٶ�
delay=2*Range_target./c0;   %�ز�ʱ��
fd1=2*Velocity_target*fc1./c0; %��Ƶ������Ƶ��
fd2=2*Velocity_target*fc2./c0; %��Ƶ������Ƶ��
%% ����ʱ����Ͷ�������
r=ones(carriers,symbols);%512�����64�����ŵľ�������������֮����ŵ�����
k_r1=ones(carriers,1);
k_r2=ones(carriers,1);
k_v1=ones(1,symbols);
k_v2=ones(1,symbols);
 %��Ƶ�µ�
for k=1:carriers
    k_r1(k,1)= exp(-1i*2*pi*(k-1)*delta_f1*delay);
end
for k=1:symbols
    k_v1(1,k)= exp(1i*2*pi*(k-1)*fd1*T1); 
end
r1=r.*(k_r1*k_v1);%��Ƶ�µ��ŵ����󣨻�û�н��п鴦��
 %��Ƶ�µ�
for k=1:carriers
    k_r2(k,1)= exp(-1i*2*pi*(k-1)*delta_f2*delay);
end
for k=1:symbols
    k_v2(1,k)= exp(1i*2*pi*(k-1)*fd2*T2); 
end
r2=r.*(k_r2*k_v2);%��Ƶ�µ��ŵ����󣨻�û�н�����״���� 
%% ��Ƶ��״����Ƶ��״ ��Ƶ�����ΪK=4��
K=4;
%% ���ڶ������ŵ����������״�Ϳ�״����
  % ���������µ��ŵ�����R3��R4
  R3=zeros(carriers,symbols);% ��Ƶ��״
  R4=zeros(carriers,symbols);% ��Ƶ��״
  R5=zeros(carriers,symbols);% ��Ƶ��״
  % ��һ������Ƶ��״����
  for i=1:K:symbols
      R3(:,i)=r1(:,i);
  end
  % �ڶ�������Ƶ��״����
  for i=1:K:carriers
      R4(i,:)=r2(i,:);
  end
  % ����������Ƶ��״����
  for i=1:K:carriers
      R5(i,:)=r1(i,:);
  end
  RMSE_V1=zeros(1,41);
  RMSE_Vidft=zeros(1,41); 
  M=500;%�ۼӴ���
  for SNR=-30:10
      RMSE1=0;
      RMSE2=0;
      tic
      for j=1:M
 %% AWGN�ŵ�
 R3_noise=awgn(R3,SNR);
 R4_noise=awgn(R4,SNR);
 R5_noise=awgn(R5,SNR);
 %% �ٶȾۺϺܼ򵥣����� exp(1i*2*pi*(k-1)*fd2*T2);���Կ�����fd2*T2=fd1*T1ʱ���ߵ�Ƶ�Ϳ���ֱ����ӣ��Ӷ��������������
 %                   ͬʱ����ѹ����֪�Կ�״�����Ż����ó���Ч���ȽϺá�ʵ�����ͨ����Ӳ�ͬCP�õ�����֪���ɲ�������Ӳ�ͬCP����ѯ�ʣ�
     %% ��һ�������Ƶ��״���ٶ�
    velocity_D1=zeros(1,symbols);
    search_V_D1=zeros(1,symbols);
    for i=1:K:carriers
        velocity_D1=R5_noise(i,:);
        D1_FFT=fft(velocity_D1);
        search_V_D1=search_V_D1+abs(D1_FFT);
    end
    search_V_D1=search_V_D1/128;
    [max_V_D1,index_V_D1]=max(search_V_D1);
    V1=(index_V_D1-1)*c0/(2*symbols*fc1*T1);
    RMSE1=RMSE1+(Velocity_target-V1)^2;
    %% �ڶ����� ��ۺϺ���ٶ�
    search_V_zuizhong=zeros(1,symbols);
    for i=1:K:512
    x1=R4_noise(i,:);
    V_fft_zuizhong=fft(x1);
    search_V_zuizhong=search_V_zuizhong+abs(V_fft_zuizhong);
    end
    G=fft(eye(symbols,symbols))/sqrt(symbols);%������ɢ�渵��Ҷ����
    G=inv(G);
    % ����G1������129:640��
    G1=zeros(symbols,symbols);
    for i=1:K:64
      G1(i,:)=G(i,:);
    end
    for i=1:512
        y=R3_noise(i,:);
        velocity_CS_D1=chat_FISTA(G1,y.',0.4);
        D1_CS_FFT=velocity_CS_D1';
        search_V_zuizhong=search_V_zuizhong+abs(D1_CS_FFT);
    end
    search_V_zuizhong=search_V_zuizhong/640;
    [max_V_zuizhong,index_V_zuizhong]=max(search_V_zuizhong(1,1:16 ));
    V_dft=(index_V_zuizhong-1)*c0/(2*symbols*fc1*T1);
    RMSE2=RMSE2+(Velocity_target-V_dft)^2;
      end
      toc
     RMSE_V1(1,SNR+31)=sqrt(RMSE1/M);
     RMSE_Vidft(1,SNR+31)=sqrt(RMSE2/M);
  end
  %% ����RMSE�ٶ�ͼ
    figure(1)%�ٶ�
    SNR_V=-30:1:10;
    semilogy(SNR_V,RMSE_V1,'b'); 
    hold on
    semilogy(SNR_V,RMSE_Vidft,'r');
    title('velocity Measurement Performance');
    xlabel('SNR/dB');
    ylabel('RMSE(v)/m/s');
    legend('High frequency comb pilot signal','CA-based staggered pilot signal')