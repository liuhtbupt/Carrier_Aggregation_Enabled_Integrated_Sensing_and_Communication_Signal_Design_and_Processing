 %% 代码说明——基于频段间载波聚合的交错导频的距离 RMSE
%代码功能：基于载波聚合的交错导频信号+改进型sturm算法进行测距测速
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
  % 第一个：高频块状处理
  for i=1:K:symbols
      R3(:,i)=r1(:,i);
  end
  % 第二个：低频梳状处理
  for i=1:K:carriers
      R4(i,:)=r2(i,:);
  end
    %% 用压缩感知
    F=ifft(eye(carriers,carriers))/sqrt(carriers);%创建离散逆傅里叶矩阵
    F=inv(F);
    F1=zeros(128,carriers);
    for i=1:128
    F1(i,:)=F(i,:);
    end  
   RMSE_R1=zeros(1,21);
   RMSE_Ridft=zeros(1,6); 
  M=500;%累加次数
  for SNR=-30:-10
      RMSE1=0;
       RMSE2=0;
      tic
      for j=1:M
      
 %% AWGN信道
 R3_noise=awgn(R3,SNR);
 R4_noise=awgn(R4,SNR);
 %% 进行聚合操作。思路：1）先将R4_noise的导频行取出，整合为128行，64列的矩阵R4_noise_zhenghe
   %                     2）对R4_noise_zhenghe利用512点的压缩感知方法求出频谱图，并与高频的频谱图相加，看是否可以提高抗噪声性能
     %% 第一步：整合
     R4_niose_zhenghe=zeros(carriers,symbols);
     m=1;
     for i=1:K:carriers 
         R4_niose_zhenghe(m,:)=R4_noise(i,:);
         m=m+1;
     end
    %% 第二步： 计算R1
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
    %% 第五步：低频整合后用512点的IDFT求距离频谱图
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
  %% 绘制RMSE距离图
    figure(1)%距离
    SNR_R=-30:1:-10;
    semilogy(SNR_R,RMSE_R1,'b'); 
    hold on
    semilogy(SNR_R,RMSE_benwen_R,'r');
    title('Range Measurement Performance');
    xlabel('SNR/dB');
    ylabel('RMSE(R)/m');
    legend('High frequency block pilot signal','CA-based staggered pilot signal')