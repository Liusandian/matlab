clear all;
close all;
j=sqrt(-1);
Array_N=4;                           % 阵元数
Array_Len=0.08;%Array_Len=0.1;
Source_N=1;                          % 信号数
Degrad=pi/180; 
Fre_Jiange=64000/1024;%64000/512; %16000/512;                  %频率间隔=62.5HZ,fs=64Khz
%%  每个频率点的间隔即频率的分辨率Ferr则是由采样率SRate和时域上所取的采样点数N来决定的,Ferr=SRate/N
%%采样率必须大于信号最高频率的2倍.现在的声卡一般采用44KHZ的采样率,44KHZ>20KHZ*2,其实就是这个意思,这里取fs=64000hz,足够大了
c=340;     %c=343;                            %声速
OverlapLen=512;%256;                 %连续段的重叠长度  取多大合适
FrameLen=1024;%512;                  %窗长
DataLen=1024;%512;                    %数据长度
Fft_Num=1024;%512;                    

Rxx=zeros(Array_N,Array_N);                 %定义协方差矩阵
Start_Fre_Point=2;                     %起始频率点，512个点傅立叶变换，起始频率约为100HZ
End_Fre_Point=40;                      %结束频率点，512个点傅立叶变换。结束频率约为3750HZ
Fre_Step=1;                          %步长
%input = load('E:\MUSIC_xianzhen512point\zhoulan_data4element_10cm_2_0m_60_1.txt');
 input = load('D:\DTV data\Hao en\single\danshenyuan\64k_five_2m_SI_Combine_20c_2.txt');
input2=input(1:400000,:); 
% input = load('E:\data_4element_xianzhen\64k\zhoulan_data4element_10cm_2_0m_30_1.txt');%走廊信号不行
 Frame_Num =length(input2)/(FrameLen-OverlapLen);%窗长减去重叠段，是否要减去1？？
Frame_Num=floor(Frame_Num);
%figure(1);
channel0=input2(:,1);
channel1=input2(:,2);
channel2=input2(:,3);
channel3=input2(:,4);
channel4=input2(:,5);
% 
figure(31);
plot(channel3);
% 
% figure(41);
i=0;
Start_Frame_Number=720;%977;%710;%1152;%42; %600(success)                  %起始帧  (500-600是合理区间)
End_Frame_Number=Frame_Num-1;%3906;%1210;%1502;%75;%%781-1
for FrameIndx=Start_Frame_Number:End_Frame_Number     %帧序数%%720-781，channel1_fftData=1×185344=181×1024
i=i+1;%init i=0;
Buf_channel1=channel1(((FrameIndx-1)*FrameLen+1-OverlapLen*(FrameIndx-1)):FrameLen*FrameIndx-OverlapLen*(FrameIndx-1))';    % 半重叠窗
Buf_channel2=channel2(((FrameIndx-1)*FrameLen+1-OverlapLen*(FrameIndx-1)):FrameLen*FrameIndx-OverlapLen*(FrameIndx-1))';
Buf_channel3=channel3(((FrameIndx-1)*FrameLen+1-OverlapLen*(FrameIndx-1)):FrameLen*FrameIndx-OverlapLen*(FrameIndx-1))';
Buf_channel4=channel4(((FrameIndx-1)*FrameLen+1-OverlapLen*(FrameIndx-1)):FrameLen*FrameIndx-OverlapLen*(FrameIndx-1))';

channel1_fftData(((i-1)*Fft_Num+1):i*Fft_Num)=fft(Buf_channel1,Fft_Num );              % 对帧数据进行1024点傅立叶变换，获得变换后的数据fftData    
channel2_fftData(((i-1)*Fft_Num+1):i*Fft_Num)=fft(Buf_channel2,Fft_Num );              %舍去第一列数据，B&K麦克风数据     
channel3_fftData(((i-1)*Fft_Num+1):i*Fft_Num)=fft(Buf_channel3,Fft_Num );                   
channel4_fftData(((i-1)*Fft_Num+1):i*Fft_Num)=fft(Buf_channel4,Fft_Num );                
    
% X(1,(((i-1)*FrameLen+1):i*FrameLen))=channel1_fftData(((i-1)*FrameLen+1):i*FrameLen);            %阵列输出信号

end
X(1,:)=channel1_fftData;
X(2,:)=channel2_fftData;
X(3,:)=channel3_fftData;
X(4,:)=channel4_fftData;
%PmusicSum=zeros(180,180);  %PmusicSum=zeros(180,180);                %定义空间谱求和矩阵
PmusicSum=zeros(1,181);  %PmusicSum=zeros(180,180);                %定义空间谱求和矩阵
Fig_Num=1;
%%for FrameIndx=2:5%Start_Frame_Number:Start_Frame_Number+5
for FrameIndx=5:7
%=========================================================================%
 for FreIndx=Start_Fre_Point:Fre_Step:End_Fre_Point     %选取频率点  ，2-40（end），Fre_step=1;       
pinlv=(FreIndx-1)*Fre_Jiange;%%62.5～2437.5Hz ,Fre_Jiange=64000/1024，0，62.5，125，...，2500Hz；（2500=62.5×40）
Y1=X(1:4,((FrameIndx-1)*Fft_Num+FreIndx));     %%FreIndx=40,FrameIndx=7,对应Y1=X(1:4,6184)     %选择对应频率点
%6184=(7-1)*1024+40;
Y2=X(1:4,(FrameIndx*Fft_Num+FreIndx));%FrameIndx=7（此时）
Y3=X(1:4,((FrameIndx+1)*Fft_Num+FreIndx));
Y4=X(1:4,((FrameIndx+2)*Fft_Num+FreIndx));
Y5=X(1:4,((FrameIndx+3)*Fft_Num+FreIndx));

Rxx1=Y1*Y1';                                         
Rxx2=Y2*Y2';  
Rxx3=Y3*Y3';
Rxx4=Y4*Y4';
Rxx5=Y5*Y5';
%*****************************第二部分：特征值分解**************************%
Rxx=(Rxx1+Rxx2+Rxx3+Rxx4+Rxx5)./5;    %输出信号相关矩阵
[U,s,v]=svd(Rxx);            %相关矩阵特征值分解
%*************************************************************************%
Vs=U(:,1:Source_N);            %信源子空间
Vn=U(:,Source_N+1:Array_N);      %噪声子空间，Array_N=4
%******************************第三部分：求空间谱****************************% 
for pitch=0%for pitch=1:180
for azimuth= 0:90%=1:180
% AA=[1;exp(-j*2*pi*pinlv*Array_Len*sin(azimuth*Degrad)*cos(pitch*Degrad)/c);
  AA=[1;exp(-j*2*pi*pinlv*Array_Len*sin(azimuth*Degrad)/c);
       exp(-j*2*pi*pinlv*3*Array_Len*sin(azimuth*Degrad)/c);
       exp(-j*2*pi*pinlv*4*Array_Len*sin(azimuth*Degrad)/c);
    ];
%   AA=[1;exp(-j*2*pi*pinlv*Array_Len*sin(azimuth*Degrad)/c);
%        exp(-j*2*pi*pinlv*2*Array_Len*sin(azimuth*Degrad)/c);
WW=AA'*Vn*Vn'*AA;
%Pmusic(pitch,azimuth)=abs(8./WW);                         %空间谱
 Pmusic(1,azimuth+91)=abs(8./WW); %为何分子为8？
 %Pmusic(1,azimuth)=abs(8./WW); 
 %Pmusic(1,azimuth)=abs(8./WW);
% Pmusic(1,pitch)=abs(8./WW);  
end %for azimuth=1:180
end %for pitch=0%pitch=1:180
 PmusicSum=PmusicSum+Pmusic;                                %空间谱求和         
 end% for FreIndx=Start_Fre_Point:Fre_Step:End_Fre_Point     %选取频率点 

 figure(Fig_Num)
 Fig_Num=Fig_Num+1;
x=-90:1:90;%x=1:1:180;
%  y=1:1:180;
plot(x,PmusicSum(1,:));
title('MUSIC空间谱');
xlabel('方位角'); 
% ylabel('俯仰角');
% zlabel('角谱');
grid on; 
i=0;

end%for FrameIndx=2:2+5%for FrameIndx=Start_Frame_Number:Start_Frame_Number+5
