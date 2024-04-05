 %% ����˵����������Ƶ�μ��ز��ۺϵĽ���Ƶ�ľ��� RMSE
%���빦�ܣ������ز��ۺϵĽ���Ƶ�ź�+�Ľ���sturm�㷨���в�����
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
  % ��һ������Ƶ��״����
  for i=1:K:symbols
      R3(:,i)=r1(:,i);
  end
  % �ڶ�������Ƶ��״����
  for i=1:K:carriers
      R4(i,:)=r2(i,:);
  end
    %% ��ѹ����֪
    F=ifft(eye(carriers,carriers))/sqrt(carriers);%������ɢ�渵��Ҷ����
    F=inv(F);
    F1=zeros(128,carriers);
    for i=1:128
    F1(i,:)=F(i,:);
    end  
   RMSE_R1=zeros(1,21);
   RMSE_Ridft=zeros(1,6); 
  M=500;%�ۼӴ���
  for SNR=-30:-10
      RMSE1=0;
       RMSE2=0;
      tic
      for j=1:M
      
 %% AWGN�ŵ�
 R3_noise=awgn(R3,SNR);
 R4_noise=awgn(R4,SNR);
 %% ���оۺϲ�����˼·��1���Ƚ�R4_noise�ĵ�Ƶ��ȡ��������Ϊ128�У�64�еľ���R4_noise_zhenghe
   %                     2����R4_noise_zhenghe����512���ѹ����֪�������Ƶ��ͼ�������Ƶ��Ƶ��ͼ��ӣ����Ƿ������߿���������
     %% ��һ��������
     R4_niose_zhenghe=zeros(carriers,symbols);
     m=1;
     for i=1:K:carriers 
         R4_niose_zhenghe(m,:)=R4_noise(i,:);
         m=m+1;
     end
    %% �ڶ����� ����R1
    search_R=zeros(carriers,1);
    for i=1:K:symbols
        x=R3_noise(:,i);
        R_ifft=ifft(x,carriers);
        search_R=search_R+abs(R_ifft);
    end
    search_R=search_R/(symbols/K);
    [max_R,index_R]=max(search_R);
    R_1=(index_R-1)*c0/(2*carriers*delta_f1);
    RMSE1=RMSE1+(Range_target-R_1)^2;
    %% ���岽����Ƶ���Ϻ���512���IDFT�����Ƶ��ͼ
     search_R_zuizhong=zeros(carriers,1);
     for i=1:symbols 
        x1=R4_niose_zhenghe(1:128,i);
        R_ifft_zuizhong=chat_FISTA(F1,x1,400);
        search_R_zuizhong=search_R_zuizhong+abs(R_ifft_zuizhong);
    end
    for i=1:K:symbols
        x2=R3_noise(:,i);
        R1_ifft=ifft(x2);
        search_R_zuizhong=search_R_zuizhong+abs(R1_ifft);
    end
    search_R_zuizhong=search_R_zuizhong/80;
    [max_R_zuizhong,index_R_zuizhong]=max(search_R_zuizhong);
    R_idft=(index_R_zuizhong-1)*c0/(2*carriers*delta_f1);
    RMSE2=RMSE2+(Range_target-R_idft)^2;
      end
      toc
      RMSE_R1(1,SNR+31)=sqrt(RMSE1/M);
      RMSE_Ridft(1,SNR+31)=sqrt(RMSE2/M);
  end
  %% ����RMSE����ͼ
    figure(1)%����
    SNR_R=-30:1:-10;
    semilogy(SNR_R,RMSE_R1,'b'); 
    hold on
    semilogy(SNR_R,RMSE_benwen_R,'r');
    title('Range Measurement Performance');
    xlabel('SNR/dB');
    ylabel('RMSE(R)/m');
    legend('High frequency block pilot signal','CA-based staggered pilot signal')