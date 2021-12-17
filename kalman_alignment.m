%-----------kalman初始对准-------------% 王新龙 惯性导航基础第二版 P254
clear 
clc
close all
%glvs

N=1000;  %迭代次数
wie=7.2916*10^-5; % 地球自传角速度
g=9.84;               %重力速度
deg=pi/180;% 1°与弧度的转换
L=45*deg;           %纬度
w=[0 wie*cos(L) wie*sin(L)];% 东北天方向上的角速度分量

% sins1 = load('D:\2021-7-19-16：21：43 RawData_1.txt');
% avp0 = avpset([0;0;0],[0;0;0],[24.618804*glv.deg;118.039968*glv.deg;0]); %设置初始导航坐标系下的avp，单位：a:deg; v:m/s; p:deg & m
% sins1_imu = [sins1(:,3:5) sins1(:,6:8)];
%  
% % 粗对准
% [attsb1, qnb1, Cnb1, eb1, db1] = alignsb(sins1_imu(2000:5000,:),avp0(7:9));  %粗对准后获得当前时刻的欧拉角、四元数、初始状态矩阵、陀螺漂移、加表零偏
% Cb1b = a2mat(attsb1(1:3));
% 
% C=Cb1b;
C=a2mat([deg deg deg]);  %假设已完成初对准

%X=[dve dvn fe fn fu dtx,dty kx ky kz]'; dve，dvn 东北-向速度误差，fe,fn, fu 失准角
%dtx,dty 加速度计零偏  kx,ky,kz 陀螺零偏
% dX=FX+BW
% Z=HX+DV
f=[0 2*w(3) 0 -g 0       
   -2*w(3) 0 g 0 0
   0 0 0 w(3) -w(2)
   0 0 -w(3) 0 0
   0 0 w(2) 0 0];
TT=[C(1,1) C(1,2) 0 0 0      
   C(2,1) C(2,2) 0 0 0
   0 0 C(1,1) C(1,2) C(1,3)
   0 0 C(2,1) C(2,2) C(2,3)
   0 0 C(3,1) C(3,2) C(3,3)];

F=[f,TT;zeros(5),zeros(5)];  %状态系数转移矩阵
%W=(zeros(10,1))';  
%-----------Kalman计算初始值设置----------%
xr=[0.1 0.1 deg deg deg zeros(1,5)]';  % 状态的真实值--初值
xe=zeros(1,10)';                       % 状态估计初值
xout=zeros(10,N);   
Q=diag([(0.00005*g)^2 (0.00005*g)^2 (0.01*deg/3600)^2 (0.01*deg/3600)^2 (0.01*deg/3600)^2 0 0 0 0 0]);% 状态噪声方差（课本上应该是零偏的平方（0.01*deg/3600)^2），应该是印刷错误）
P=diag([0.1^2 0.1^2 deg^2 deg^2 deg^2 (0.0001*g)^2 (0.0001*g)^2 (0.02*deg/3600)^2 (0.02*deg/3600)^2 (0.02*deg/3600)^2]); %初始协方差矩阵
R=diag([0.1^2 0.1^2]);                 % 量测噪声方差阵
H=[1 0 0 0 0 0 0 0 0 0    
   0 1 0 0 0 0 0 0 0 0];               % 量测转移矩阵
xrk_1=xr;
xk_1=xe;
pk_1=P;

%------- Kalman离散化------------------%
A=F;
B=eye(10);
D=zeros(2,10);                         % 这里B和D实际上都没有用，只是为了能够使用ss函数才设成D=zeros(2,10)
G=ss(A,B,H,D);
T=1;
Gd=c2d(G,T);                           % 连续模型转离散化的函数
FF=Gd.a;                               % 离散化的转移矩阵



%----------Kalman计算过程-----------%
for i =1:N
    
    xrk=A*xrk_1+(mvnrnd(zeros(10,1),Q,1))';%这里计算实际的状态究竟是用离散化之后的模型还是用连续模型呢，用连续模型得到的效果还可以，离散的不行？
    z=H*xrk_1+(mvnrnd([0;0],R,1))';
    %一步预测
    xk_=FF*xk_1;
    pk_=FF*pk_1*FF'+Q;
    %更新
    K=pk_*H'/(H*pk_*H'+R); %卡尔曼增益
    xk=xk_+K*(z-H*xk_);    %估计值
    pk=(eye(10)-K*H)*pk_;  %协方差估计值
    xout(:,i)=xk;
    xk_1=xk;
    pk_1=pk;
    xrk_1=xrk;
end
%--------画图----------%
t=1:N;
%x1=xout'
figure (1);
plot(t,xout(1,:),t,xout(2,:));
xlabel('t/s');
ylabel('m');
title('速度误差')
figure (2)
plot(t,xout(3,:)*60*180/pi,t,xout(4,:)*60*180/pi);%化成分（角）输出
xlabel('t/s');
ylabel('arcmin');
title('东向、北向失准角');
figure (3)
plot(t,xout(5,:)*60*180/pi);
xlabel('t/s');
ylabel('arcmin');
title('天向失准角')




