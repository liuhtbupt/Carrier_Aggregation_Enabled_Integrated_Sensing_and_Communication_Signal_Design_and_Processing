%% 代码说明——基于频段间载波聚合的交错导频测速
%代码功能：
%编写人：刘浩田
%日期：2022.11.30
%% 清除工作区
clear;
clc;
%% 高频部分系统参数
fc1=24e9;	%高频载波频率
symbols=64; %OFDM符号64
carriers=512; %子载波个数
delta_f1=120e3;% 高频子载波间隔
c0=3e8;  %光速
cp = 1.33*10^(-6);%满足最大测距200米
T1=1/delta_f1+cp;%高频符号长度
%% 低频部分系统参数
fc2=5.9e9;	%低频载波频率
delta_f2=30e3; %低频子载波间隔
T2=1/delta_f2+5.975141242935*10^(-6);%低频符号长度
%% 目标参数
Range_target=117;    %目标距离
Velocity_target=30;  %目标速度
delay=2*Range_target./c0;   %回波时延
fd1=2*Velocity_target*fc1./c0; %高频多普勒频移
fd2=2*Velocity_target*fc2./c0; %低频多普勒频移
%% 构建时延项和多普勒项
r=ones(carriers,symbols);%512点采样64个符号的经过分素数除法之后的信道矩阵
k_r1=ones(carriers,1);
k_r2=ones(carriers,1);
k_v1=ones(1,symbols);
k_v2=ones(1,symbols);
 %高频下的
for k=1:carriers
    k_r1(k,1)= exp(-1i*2*pi*(k-1)*delta_f1*delay);
end
for k=1:symbols
    k_v1(1,k)= exp(1i*2*pi*(k-1)*fd1*T1); 
end
r1=r.*(k_r1*k_v1);%高频下的信道矩阵（还没有进行块处理）
 %低频下的
for k=1:carriers
    k_r2(k,1)= exp(-1i*2*pi*(k-1)*delta_f2*delay);
end
for k=1:symbols
    k_v2(1,k)= exp(1i*2*pi*(k-1)*fd2*T2); 
end
r2=r.*(k_r2*k_v2);%低频下的信道矩阵（还没有进行梳状处理） 
%% 高频块状、低频梳状 导频间隔均为K=4；
K=4;
%% 现在对两个信道矩阵进行梳状和块状处理
  % 创建两个新的信道矩阵R3，R4
  R3=zeros(carriers,symbols);% 高频块状
  R4=zeros(carriers,symbols);% 低频梳状
  R5=zeros(carriers,symbols);% 高频梳状
  % 第一个：高频块状处理
  for i=1:K:symbols
      R3(:,i)=r1(:,i);
  end
  % 第二个：低频梳状处理
  for i=1:K:carriers
      R4(i,:)=r2(i,:);
  end
  % 第三个：高频梳状处理
  for i=1:K:carriers
      R5(i,:)=r1(i,:);
  end
  RMSE_V1=zeros(1,41);
  RMSE_Vidft=zeros(1,41); 
  M=500;%累加次数
  for SNR=-30:10
      RMSE1=0;
      RMSE2=0;
      tic
      for j=1:M
 %% AWGN信道
 R3_noise=awgn(R3,SNR);
 R4_noise=awgn(R4,SNR);
 R5_noise=awgn(R5,SNR);
 %% 速度聚合很简单，根据 exp(1i*2*pi*(k-1)*fd2*T2);可以看到当fd2*T2=fd1*T1时，高低频就可以直接相加，从而可以提升信噪比
 %                   同时再用压缩感知对块状进行优化，得出的效果比较好。实现相等通过添加不同CP得到，不知道可不可以添加不同CP（待询问）
     %% 第一步：求高频梳状的速度
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
    %% 第二步： 求聚合后的速度
    search_V_zuizhong=zeros(1,symbols);
    for i=1:K:512
    x1=R4_noise(i,:);
    V_fft_zuizhong=fft(x1);
    search_V_zuizhong=search_V_zuizhong+abs(V_fft_zuizhong);
    end
    G=fft(eye(symbols,symbols))/sqrt(symbols);%创建离散逆傅里叶矩阵
    G=inv(G);
    % 构造G1，用于129:640行
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
  %% 绘制RMSE速度图
    figure(1)%速度
    SNR_V=-30:1:10;
    semilogy(SNR_V,RMSE_V1,'b'); 
    hold on
    semilogy(SNR_V,RMSE_Vidft,'r');
    title('velocity Measurement Performance');
    xlabel('SNR/dB');
    ylabel('RMSE(v)/m/s');
    legend('High frequency comb pilot signal','CA-based staggered pilot signal')