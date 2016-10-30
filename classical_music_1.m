%MUSIC ALOGRITHM
%DOA ESTIMATION BY CLASSICAL_MUSIC
clear all;
%close all;
clc;

source_number=3;%信元数
sensor_number=8;%阵元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=[pi/4 pi/6].';%信号频率
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;%信号波长  
d=0.5*l;%阵元间距
snr=0;%信噪比

source_doa=[-45 60];%两个信号的入射角度
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l)].';%阵列流型

s=sqrt(10.^(snr/10))*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号
powersignal=[10 10];
for i=1:2
    for ii=1:N_x
        S(i,ii)=sqrt(0.5*10^(powersignal(i)/10))+j*sqrt(0.5*10^(powersignal(i)/10));
    end
end
x=A*S;
xx=x+1/sqrt(2)*(randn(8,N_x)+j*randn(8,N_x));



R=x*x'/snapshot_number;
disp(R);
%[V,D]=eig(R);
%Un=V(:,1:sensor_number-source_number);
%Gn=Un*Un';
[U,S,V]=svd(R);
disp(S);
Un=U(:,source_number+1:sensor_number);
Gn=Un*Un';

searching_doa=-90:0.1:90;%线阵的搜索范围为-90~90度
 for i=1:length(searching_doa)
   a_theta=exp(-j*(0:sensor_number-1)'*2*pi*d*sin(pi*searching_doa(i)/180)/l);
   Pmusic(i)=1./abs((a_theta)'*Gn*a_theta);
 end
plot(searching_doa,10*log10(Pmusic),'r');
%axis([-90 90 -90 90]);
xlabel('入射角/degree');
ylabel('空间谱/dB');
legend('Music Spectrum');
title('经典MUSIC估计');
grid on;
