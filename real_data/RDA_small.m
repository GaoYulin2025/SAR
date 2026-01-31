%%卫星轨道参数设置
H = 755e3; %卫星轨道高度
% 卫星轨道速度Vr计算
EarthMass = 6e24; %地球质量(kg)
EarthRadius = 6.37e6; %地球半径6371km
Gravitational = 6.67e-11; %万有引力常量
% 姿态参数
phi = 20 * pi / 180; % 俯仰角+20°
incidence = 20.5 * pi / 180; % 入射角
%计算等效雷达速度(卫星做圆周运动的线速度)
Vr = sqrt(Gravitational*EarthMass/(EarthRadius + H)); %第一宇宙速度

%景中心斜距R_eta_c和最近斜距R0
% 斜视角theta_rc由轨道高度、俯仰角、入射角计算得出
R_eta_c = H / cos(incidence); %景中心斜距
R0 = H / cos(phi);
theta_r_c = acos(R0/R_eta_c); %斜视角，单位为弧度;斜视角为4.6°

%%信号参数设置
%   电磁波参数
c = 3e+8; % 光速
Vs = Vr; % 卫星平台速度
Vg = Vr; % 波束扫描速度
La = 15; %方位向天线长度->椭圆的长轴
Lr=1.5;%距离向天线尺寸——>椭圆的短轴
f0 = 5.4e+9; % 雷达工作频率
lambda = c / f0; %电磁波波长

%  距离向信号参数
Tr = 40e-6; % 发射脉冲时宽
Br = 2.8 * 6e6; % 距离向信号带宽
Kr = Br / Tr; % 距离向调频率
alpha_os_r = 1.2; % 距离过采样率
Nrg = 2500; % 距离线采样点数
Fr = alpha_os_r * Br; % 距离向采样率

%  方位向信号参数
alpha_os_a = 1.7; % 方位过采样率(高过采样率避免鬼影目标)
Naz = 1600; % 距离线数
delta_f_dop = 2 * 0.886 * Vr * (cos(theta_r_c)) / La; % 多普勒带宽
Fa = alpha_os_a * delta_f_dop; % 方位向采样率
Ta = 0.886 * lambda * R_eta_c / (La * Vg * cos(theta_r_c)); %目标照射时间

%  景中心点(原点)的参数
time_eta_c = -R_eta_c * sin(theta_r_c) / Vr; % 波束中心穿越时刻
f_eta_c = 2 * Vr * sin(theta_r_c) / lambda; % 多普勒中心频率

%  合成孔径参数
rho_r = c / (2 * Fr); % 距离向分辨率
rho_a = Vr/Fa; % 距离向分辨率|La / 2
theta_bw = 0.886 * lambda / La; % 方位向3dB波束宽度




%%时间轴参数
%设置时间采样间隔，频率采样间隔
%距离向，方位向总的时间
Trg=Nrg/Fr;
Taz=Naz/Fa;
%距离向，方位向时间采样间隔
Gap_t_tau=1/Fr;
Gap_t_eta=1/Fa;
%距离向，方位向频率采样间隔
Gap_f_tau=Fr/Nrg;
Gap_f_eta=Fa/Naz;
%生成距离向，方位向时间轴
time_tau_r=2*R_eta_c/c+(-Trg/2:Gap_t_tau:Trg/2-Gap_t_tau);%大小为Nrg
time_eta_a=time_eta_c+(-Taz/2:Gap_t_eta:Taz/2-Gap_t_eta);%大小为Naz
[Ext_time_tau_r, Ext_time_eta_a] = meshgrid(time_tau_r, time_eta_a); % 设置距离时域-方位时域二维网络坐标，其中Ext_time_tau_r,大小为Naz*nrg,按行重复，Ext_time_eta_a大小为naz*nrg按列重复
%生成距离向，方位向频率轴
f_tau=(-Fr/2:Gap_f_tau:Fr/2-Gap_f_tau);
f_eta=(-Fa/2:Gap_f_eta:Fa/2-Gap_f_eta);
% f_tau=f_tau-(round(f_tau/Fr))/Fr;%混叠方程，混叠后的
% f_eta=f_eta-(round(f_eta/Fa))/Fa;
[Ext_f_tau, Ext_f_eta] = meshgrid(f_tau, f_eta); % 设置频率时域-方位频域二维网络坐标
%生成随距离向时间变化的最短距离
R0_tau_r=(c*time_tau_r/2)*cos(theta_r_c);
Ext_R0_tau_r=repmat(R0_tau_r,Naz,1);%扩展成Naz*Nrg的大小，其中每一行都是第一行的重复，原因是最近斜距只和距离向时间有关
%景中心的绝对频率
f_eta_c_abs = 2 * Vr * sin(theta_r_c) / lambda;
%计算模糊数
N_amb = round(f_eta_c_abs / Fa);
%计算绝对方位向频率
f_eta_abs = f_eta + N_amb * Fa;
[Ext_f_tau, Ext_f_eta_abs] = meshgrid(f_tau, f_eta_abs);

%%点目标设置
%生成三个点目标,A C同一波束中心穿越时刻，B C方位向下相同
xa=0;ya=0;
xb=ya+500;yb=xb+500;
xc = H*tan(phi+theta_bw/2)-R0*sin(phi);%计算的C点距离向坐标
yc=ya+500;
Position_x_r=[xa,xb,xc];
Position_y_a=[ya,yb,yc];

%%回波生成
Target_num = 3; %目标数量
S_echo = zeros(Naz, Nrg);
for i=1:Target_num
    %分别生成最近斜距，波束中心穿越时刻，瞬时斜距，距离向包络，方位向包络，相位信息

    %最近斜距
    R0_Target=sqrt(H^2+(R_eta_c*sin(phi)+Position_x_r(i))^2);
    %波束中心穿越时刻
    time_eta_c_Target=(Position_y_a(i)-R0_Target*tan(theta_r_c))/Vr;
    %瞬时斜距,算出了全部方位向时间的瞬时斜距
    R_eta=sqrt((R0_Target^2)+(Vr^2)*(Ext_time_eta_a-Position_y_a(i)/Vr).^2);
    %距离向包络
    Wr=(abs(Ext_time_tau_r-2*R_eta/c)<=Tr/2);
    %方位向包络
    Wa=(abs(Ext_time_eta_a-time_eta_c_Target)<=Ta/2);
    %相位信息
    Phase = exp(-1j*4*pi*f0*R_eta/c) .* exp(+1j*pi*Kr*(Ext_time_tau_r - 2 * R_eta / c).^2);
    %拼接接收的信号
    S_echo_Target=Wr.*Wa.*Phase;

    %对所有的接收信号进行叠加
    S_echo=S_echo + S_echo_Target;

end

%%时间轴矫正
S_echo = S_echo .* exp(-1j*2*pi*f_eta_c*Ext_time_eta_a); %多普勒中心校正

%%距离压缩
%定义频域滤波器
Hf=(abs(Ext_f_tau)<=Br/2).*exp(+1j*pi*Ext_f_tau.^2/Kr);
win_rng = hann(Nrg).';                 % Hann 窗
Win2D = repmat(win_rng, Naz, 1);
Hf = Hf .* Win2D;    
%距离向傅里叶变换
s1_ftau_eta=fftshift(fft(fftshift(S_echo,2),Nrg,2),2);
%匹配滤波
s1_ftau_eta=s1_ftau_eta.*Hf;
%距离向傅里叶逆变换
s1_tau_eta=fftshift(ifft(fftshift(s1_ftau_eta,2),Nrg,2),2);


%%距离徙动矫正
%距离徙动矫正在距离多普勒域，首先进行方位向傅里叶变换，然后定义徙动量
%方位向傅里叶变换
s2_tau_feta=fftshift(fft(fftshift(s1_tau_eta,1),Naz,1),1);
%定义徙动量,这里的最短斜距R0是景中心对应的最短斜距
delta_R=(((lambda * Ext_f_eta_abs).^2) .* R0) ./ (8 * (Vr^2));
%这里采用相位乘法器进行RCMC，这种方法在二维频域实现，首先定义相位乘法器
G_rcmc=exp(+1j*4*pi*Ext_f_tau.*delta_R./c);
%进行距离向傅里叶变换，转换到二维频域
s3_ftau_feta=fftshift(fft(fftshift(s2_tau_feta,2),Nrg,2),2);
%信号与相位乘法器相乘
s3_ftau_feta=s3_ftau_feta.*G_rcmc;


%距离向傅里叶逆变换
s3_tau_feta_RCMC = fftshift(ifft(fftshift(s3_ftau_feta, 2), Nrg, 2), 2); %距离向傅里叶逆变换

%%方位压缩
%首先计算方位向调频率，然后根据调频率涉及方位向匹配滤波器
%方位向调频率是根据最近斜距不断变化的
ka = 2 * Vr^2 * (cos(theta_r_c)^2) ./ (lambda .* Ext_R0_tau_r);
%方位向匹配滤波器
W_azimuth = hamming(Naz);
Haz=exp(-1j*pi.*Ext_f_eta.^2./ka).*W_azimuth;
%Offset = exp(-1j*2*pi*Ext_f_eta.*time_eta_c);%偏移滤波器，将原点搬移到Naz/2的位置，校准坐标
%回波信号与滤波器相乘
%s4_tau_feta=s3_tau_feta_RCMC.*Haz.*Offset;
s4_tau_feta=s3_tau_feta_RCMC.*Haz;
%转回二维时域
S_Image = fftshift(ifft(fftshift(s4_tau_feta, 1), Naz, 1), 1);









%%绘图


figure('name', "成像处理流程", 'Position', [100, 100, 1200, 800])

subplot(2, 2, 1);
imagesc(real(S_echo)); 
title('1. 原始回波 (实部)');
xlabel('距离向'); ylabel('方位向');

subplot(2, 2, 2);
imagesc(abs(s1_tau_eta)); 
title('2. 距离压缩后 (幅度)');
xlabel('距离向'); ylabel('方位向');

subplot(2, 2, 3);
imagesc(abs(s3_tau_feta_RCMC)); 
title('3. RCMC校正后 (距离多普勒域幅度)');
xlabel('距离向'); ylabel('方位向频率');

subplot(2, 2, 4);
imagesc(abs(S_Image)); 
title('4. 最终方位压缩结果 (幅度)');
xlabel('距离向'); ylabel('方位向');
colormap('jet'); % 使用彩色图更清晰





