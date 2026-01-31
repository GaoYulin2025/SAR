%% RDA_High_Distinction.m
% 核心：通过 Gamma 校正和 S 型截断，最大化黑白区分度
% 改进版: 使用标准的 load_AGC_block / load_DATA_block 读取数据
clear; clc; close all;

%% 1. 加载与核心处理
fprintf('读取数据...\n');

% 1.1 加载参数
if exist('CD_run_params.mat', 'file')
    load('CD_run_params.mat');
else
    error('缺参数文件: CD_run_params.mat'); 
end

% 1.2 设置读取块参数 (参考 run_RD_simulation.m)
Block_Index = 1;      
fprintf('正在加载块 %d 数据...\n', Block_Index);

% 1.3 使用辅助函数读取数据
% 注意: 必须确保 load_AGC_block 和 load_DATA_block 在路径中 (project 目录下)
% 第一个参数传空字符串 ''，因为 UseMATfiles=1 时不需要前缀，或者由 params 控制
AGC_values = load_AGC_block('', first_rg_line, Nrg_lines_blk, Block_Index, UseMATfiles);
data = load_DATA_block('', output_path, Nrg_lines_blk, Nrg_cells, AGC_values, Block_Index, UseMATfiles);

[Na_orig, Nr_orig] = size(data);    %获取数据本身的尺寸   Na是方位向  慢时间  Nr是距离向  快时间
data = double(data);

C = c;                  %光速
Lambda = C / f0;        %波长 
Fs = Fr;                %距离向采样频率
Dur = Tr;               %时宽  雷达发射时间
PRFo = PRF;             %脉波频率   方位向采样频率
dt = 1/Fs;              %距离向采样时间间隔
fast_t = (0:Nr_orig-1)*dt + 2*R0/C;     %快事件轴   距离时间轴
Ra = fast_t * C / 2;    %斜距   目标和雷达的距离
V_r = 7062;             %有效雷达速度
Ambiguity_Num = -6;     %多普勒模糊数
C0 = 655; C1 = -14.7; C2 = 0.05;     %多普勒中心频率的积分式系数   这个是给那个你分完块之后的那个坐标系做准备的

%% 2. 距离压缩 (稳妥修正版)
fprintf('1/4 距离压缩...\n');

% 1. 计算参数
N_ref = floor(Dur * Fs);      %floor函数作用是取整    距离向的采样点数  有多少个  
% 确保 FFT 长度足够长，防止混叠
Nfft_r = 2^nextpow2(Nr_orig + N_ref);    %nextpow2函数作用是  取2的幂次方

% 2. 构造参考信号 (Chirp)
t_ref = linspace(-Dur/2, Dur/2, N_ref);
Ht_r = exp(1j * pi * (-abs(Kr)) * t_ref.^2);

% --- 修正点：在时域正确加窗 ---
% 目的：只对参考信号本身加窗，而不是对补零后的频谱加窗
% 这样可以保证主瓣能量不损失，同时压低旁瓣
win_t = ones(N_ref,1).';     %这个是一行 N_ref列     .'是转置
Ht_r_weighted = Ht_r .* win_t;

% 3. 生成频域匹配滤波器
% 注意：这里不加 fftshift，保持和你原代码一致的相位中心，防止目标移位
Hw_r = conj(fft(Ht_r_weighted, Nfft_r));    %conj取复共轭  将每一行补到Nfft_r这么长

S_rg = zeros(Na_orig, Nr_orig);
for i = 1:Na_orig
    % A. 数据 FFT
    Sig_f = fft(data(i, :), Nfft_r);
    
    % B. 匹配滤波
    Sig_compressed = ifft(Sig_f .* Hw_r);
    
    % C. 截断 (保持原代码逻辑，直接截取前段)
    % 这样目标的位置就不会乱跑了
    S_rg(i, :) = Sig_compressed(1:Nr_orig); 
end

fprintf('距离压缩完成。\n');

%% 3. RCMC
fprintf('2/4 RCMC...\n');
% 对距离压缩后的数据 S_rg，在方位向（维度1）做 FFT，并把零频移到中间。
% 此时数据进入了 "距离-多普勒域" (Range-Doppler Domain)。
% 这是 RCMC 必须进行的域，因为距离弯曲的大小取决于多普勒频率。
S_rd = fftshift(fft(S_rg, [], 1), 1);        
freq_az = linspace(-PRFo/2, PRFo/2, Na_orig).';     %生成频率轴  为后面弄绝对多普勒频率做准备  
S_rcmc = zeros(Na_orig, Nr_orig);                   %储存数据用的
center_col = round(Nr_orig/2);                      %多普勒中心频率  只计算了中间那一列
%它是连接"当前像素位置"和"参数模型"的接口。它保证了你输入的像素尺寸（0~9），与计算出$C_0, C_1, C_2$时使用的测量尺是一致的。
d_norm = center_col / Nr_orig * 9;                  %这个9是有一个文件中的  他是分成了3x3的块  然后是9 划分成了0-9的坐标系   归一化距离 它负责把图像中的像素坐标，翻译成逻辑坐标
f_c = (C2*d_norm^2 + C1*d_norm + C0) + Ambiguity_Num * PRFo;   %这个就是数据的参考多普勒中心
for i = 1:Na_orig                                        
    if mod(i, 2000)==0, drawnow; end
    f_eta = freq_az(i) + f_c;        %这个是数据的绝对多普勒频率
    val = 1 - (Lambda * f_eta / (2*V_r))^2;     
    D_f = sqrt(max(val, 0));          %这一步和上面那一步是在计算徙动因子  是一个小于一的数  直视距离与斜视距离之间的比例关系
    R_mig = Ra ./ D_f;                %我们想要一个直的图像距离Ra   但在这个弯曲图像里能量被拉伸到了R_mig
    col_idx = (R_mig - Ra(1)) / (C/2 * dt) + 1;   %把这个拉伸的距离换算成了阵列的列索引
    %interp1就像一只手。我对它说："对于第1个格子，你去原图的col_idx(1)位置取数据；对于第2个格子，你去col_idx(2)位置取数据……"
    %由于col_idx往往是小数（比如512.5），所以需要'spline'（样条插值）来准确提示出那个小数位置的值。
    S_rcmc(i, :) = interp1(1:Nr_orig, S_rd(i, :), col_idx, 'spline', 0);
end

%% 4. 方位压缩
fprintf('3/4 方位压缩...\n');
S_image = zeros(Na_orig, Nr_orig);      
win_az = ones(Na_orig,1); 
for k = 1:Nr_orig
    if mod(k, 1000)==0, fprintf('   Azimuth %.0f%%\n', k/Nr_orig*100); drawnow; end
    % 距离 R_curr 越远，Ka 越小，意味着合成孔径越长，聚焦越难。
    R_curr = Ra(k);      %这一方位上的斜距  也就是这一个方位上的目标与雷达的距离
    d_n = k / Nr_orig * 9;     %这个是算归一化距离
    % 计算当前列的 "多普勒中心频率 f_dc"
    % 这决定了信号频谱的中心位置 (Squint Angle)
    f_val = (C2*d_n^2 + C1*d_n + C0) + Ambiguity_Num * PRFo;
    Ka = 2 * V_r^2 / (Lambda * R_curr);      %方位向调频率
    H_az = exp(-1j * pi * (freq_az+f_val).^2 / Ka);    %方位向匹配滤波器
    % 加上 ifftshift 确保相位正确归位
    S_image(:, k) = ifft(S_rcmc(:, k) .* H_az .* win_az);
end

%% 5. 【核心优化】高区分度增强显示
fprintf('4/4 增强显示效果...\n');

% 1. 转置与取模
img = abs(S_image); 

% 2. 对数变换
img_log = 20 * log10(img + 1e-6);

% 3. 去噪 (中值滤波 - 保留边缘)
img_clean = medfilt2(img_log, [3 3]); 

% 4. 【直方图截断 (更严格)】
sorted_vals = sort(img_clean(:));
% 丢弃底部 8% (让海水彻底变黑)
min_val = sorted_vals(floor(length(sorted_vals) * 0.08)); 
% 丢弃顶部 0.5% (只保留极强点)
max_val = sorted_vals(floor(length(sorted_vals) * 0.995)); 

% 归一化到 0~1
img_norm = (img_clean - min_val) / (max_val - min_val);
img_norm(img_norm < 0) = 0;
img_norm(img_norm > 1) = 1;

% 5. 【Gamma 校正】(拉开暗部层次)
% Gamma > 1 会压暗中间调，Gamma < 1 会提亮中间调
% 这里用 1.3 让背景更黑，目标更突出
gamma_val = 1.3; 
img_final = img_norm .^ gamma_val;

% 6. 【S型曲线增强】(可选，进一步增加"硬"度)
% 如果觉得不够硬，取消下面两行的注释
% gain = 10; cutoff = 0.4;
% img_final = 1 ./ (1 + exp(-gain * (img_final - cutoff)));

% 绘图
figure('Name', 'RADARSAT-1 高区分度版', 'Color', 'white', 'Position', [50 50 800 900]);
imagesc(img_final); 
colormap('gray'); 
axis image; axis xy;
xlabel('方位向 (Azimuth)'); ylabel('距离向 (Range)');
title('RD算法成像');

fprintf('处理完成！现在的图片应该黑白分明。\n');

%% --- 综合对比处理：修正方向版 ---
fprintf('开始生成对比图像 (修正方向)...\n');

Img_Raw = abs(S_image); 

% =========================================================
% 方案一：4x4 多视处理
% =========================================================
[rows, cols] = size(Img_Raw);
max_r = floor(rows/4)*4;
max_c = floor(cols/4)*4;
S_trim = Img_Raw(1:max_r, 1:max_c); 

blk_vals = zeros(max_r/4, max_c/4);
for r = 1:4
    for c = 1:4
        blk_vals = blk_vals + S_trim(r:4:end, c:4:end);
    end
end
ML_Strong = blk_vals / 16; 

% 增强
clim_val = prctile(ML_Strong(:), 99);
ML_Strong(ML_Strong > clim_val) = clim_val;
ML_Strong = (ML_Strong / clim_val) .^ 0.5;

% --- 绘图 1 (修正转置) ---
figure('Name', '方案一：4x4 多视结果'); 
% 【核心修改】ML_Strong.' -> 加���转置符号，把矩阵竖起来
imagesc(ML_Strong); 
colormap gray; 
axis image; axis xy; 
xlabel('方位向 (Azimuth)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('距离向 (Range)', 'FontSize', 12, 'FontWeight', 'bold');   
title('方案一：4x4 多视 (方向已修正)');


% =========================================================
% 方案二：Lee 滤波
% =========================================================
fprintf('正在进行 Lee 滤波...\n');

Img_Original = Img_Raw;
limit = prctile(Img_Original(:), 99);
Img_Original(Img_Original > limit) = limit;
Img_Original = Img_Original / limit;

% Lee 滤波
Img_Lee = my_lee_filter(Img_Original, [5 5]); 
Img_Lee_Final = Img_Lee .^ 0.5;

% --- 绘图 2 (修正转置) ---
figure('Name', '方案二：Lee 滤波结果');
% 【核心修改】Img_Lee_Final.' -> 加了转置符号
imagesc(Img_Lee_Final);
colormap gray; 
axis image; axis xy; 
xlabel('方位向 (Azimuth)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('距离向 (Range)', 'FontSize', 12, 'FontWeight', 'bold');   
title('方案二：Lee 滤波 (方向已修正)');


% =========================================================
% Lee 滤波函数 (保持不变)
% =========================================================
function out = my_lee_filter(img, window_size)
    img = double(img);
    h = ones(window_size) / prod(window_size);
    local_mean = conv2(img, h, 'same');
    local_sqr_mean = conv2(img.^2, h, 'same');
    local_var = local_sqr_mean - local_mean.^2;
    local_var(local_var < 0) = 0;
    overall_noise_var = mean(local_var(:)); 
    k = (local_var - overall_noise_var) ./ (local_var + 1e-6); 
    k(k < 0) = 0; k(k > 1) = 1;
    out = local_mean + k .* (img - local_mean);
end