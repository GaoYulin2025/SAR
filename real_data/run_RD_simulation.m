%% RDA_Raw_Only.m
% 仅展示原始单视图像 (Single Look)，不包含多视处理和滤波
% 保留基础的 Log 和 Gamma 增强以确保图像可见
clear; clc; close all;

%% 1. 加载数据
fprintf('读取数据...\n');

% 1.1 加载参数
if exist('CD_run_params.mat', 'file')
    load('CD_run_params.mat');
else
    error('缺参数文件: CD_run_params.mat'); 
end

% 1.2 设置读取块参数
Block_Index = 1;      
fprintf('正在加载块 %d 数据...\n', Block_Index);

% 1.3 读取数据
% 确保 load_AGC_block 和 load_DATA_block 在路径中
AGC_values = load_AGC_block('', first_rg_line, Nrg_lines_blk, Block_Index, UseMATfiles);
data = load_DATA_block('', output_path, Nrg_lines_blk, Nrg_cells, AGC_values, Block_Index, UseMATfiles);

[Na_orig, Nr_orig] = size(data);    
data = double(data);

% 参数设置
C = c;                  
Lambda = C / f0;        
Fs = Fr;                
Dur = Tr;               
PRFo = PRF;             
dt = 1/Fs;              
fast_t = (0:Nr_orig-1)*dt + 2*R0/C;     
Ra = fast_t * C / 2;    
V_r = 7062;             
Ambiguity_Num = -6;     
C0 = 655; C1 = -14.7; C2 = 0.05;     

%% 2. 距离压缩
fprintf('1/3 距离压缩...\n');

N_ref = floor(Dur * Fs);      
Nfft_r = 2^nextpow2(Nr_orig + N_ref);    

t_ref = linspace(-Dur/2, Dur/2, N_ref);
Ht_r = exp(1j * pi * (-abs(Kr)) * t_ref.^2);
win_t = ones(N_ref,1).';     
Ht_r_weighted = Ht_r .* win_t;

Hw_r = conj(fft(Ht_r_weighted, Nfft_r));    

S_rg = zeros(Na_orig, Nr_orig);
for i = 1:Na_orig
    Sig_f = fft(data(i, :), Nfft_r);
    Sig_compressed = ifft(Sig_f .* Hw_r);
    S_rg(i, :) = Sig_compressed(1:Nr_orig); 
end

%% 3. RCMC
fprintf('2/3 RCMC...\n');
S_rd = fftshift(fft(S_rg, [], 1), 1);        
freq_az = linspace(-PRFo/2, PRFo/2, Na_orig).';     
S_rcmc = zeros(Na_orig, Nr_orig);                   
center_col = round(Nr_orig/2);                      
d_norm = center_col / Nr_orig * 9;                  
f_c = (C2*d_norm^2 + C1*d_norm + C0) + Ambiguity_Num * PRFo;   

for i = 1:Na_orig                                        
    f_eta = freq_az(i) + f_c;        
    val = 1 - (Lambda * f_eta / (2*V_r))^2;     
    D_f = sqrt(max(val, 0));          
    R_mig = Ra ./ D_f;                
    col_idx = (R_mig - Ra(1)) / (C/2 * dt) + 1;   
    S_rcmc(i, :) = interp1(1:Nr_orig, S_rd(i, :), col_idx, 'spline', 0);
end

%% 4. 方位压缩
fprintf('3/3 方位压缩...\n');
S_image = zeros(Na_orig, Nr_orig);      
win_az = ones(Na_orig,1); 
for k = 1:Nr_orig
    R_curr = Ra(k);      
    d_n = k / Nr_orig * 9;     
    f_val = (C2*d_n^2 + C1*d_n + C0) + Ambiguity_Num * PRFo;
    Ka = 2 * V_r^2 / (Lambda * R_curr);      
    H_az = exp(-1j * pi * (freq_az+f_val).^2 / Ka);    
    S_image(:, k) = ifft(S_rcmc(:, k) .* H_az .* win_az);
end

%% 5. 图像显示 (基础增强)
fprintf('生成最终图像...\n');

% 1. 取模
img = abs(S_image); 

% 2. 对数变换 (20*log10)
% 必须步骤：否则图像除了强点外全是黑的
img_log = 20 * log10(img + 1e-6);

% 3. 直方图截断与归一化
% 简单的自动对比度拉伸，去掉极亮和极暗的噪点
limit_min = prctile(img_log(:), 5);  % 丢弃底部 5%
limit_max = prctile(img_log(:), 99); % 丢弃顶部 1%

img_norm = (img_log - limit_min) / (limit_max - limit_min);
img_norm(img_norm < 0) = 0;
img_norm(img_norm > 1) = 1;

% 4. Gamma 校正
% 提亮暗部细节
img_final = img_norm .^ 1.2; 

% 绘图
figure('Name', '原始单视 SAR 图像', 'Color', 'white');

% 注意：这里使用了转置 (.')，让横轴显示为方位向，竖轴显示为距离向，
% 这通常更符合 SAR 瀑布图的观察习惯。如果您习惯横着看，可以去掉 .'
imagesc(img_final.'); 

colormap('gray'); 
axis image; axis xy;
xlabel('方位向 (Azimuth)'); 
ylabel('距离向 (Range)');
title('原始单视复数图像 (SLC) - 仅做 Log/Gamma 增强');

fprintf('显示完成。\n');