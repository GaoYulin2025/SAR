% run_RD_simulation.m
% 基于提取的数据进行 RD 算法成像仿真
clear; clc; close all;

% 1. 加载参数 (由 extract_data 生成)
load('CD_run_params.mat');

% 调整参数以适应仿真
% 注意：为了快速测试，这里只演示处理第 1 个块的数据
% 如果想处理全图，需要将结果拼接起来
Block_Index = 1; 

% 2. 准备路径和文件名
file_pre = fullfile(output_path, ['CDdata' num2str(Block_Index)]); % 注意：extract_data保存的文件名可能没有前缀，直接是CDdataX

% 3. 读取数据 (使用项目提供的辅助函数)
fprintf('正在加载块 %d 数据...\n', Block_Index);
% 加载 AGC 值
AGC_values = load_AGC_block('', first_rg_line, Nrg_lines_blk, Block_Index, UseMATfiles);
% 加载 I/Q 数据
raw_data_block = load_DATA_block('', output_path, Nrg_lines_blk, Nrg_cells, AGC_values, Block_Index, UseMATfiles);

% raw_data_block 现在的维度应该是 [Nrg_lines_blk x Nrg_cells]
% 行是方位向(Azimuth)，列是距离向(Range)
[Na, Nr] = size(raw_data_block);

% -----------------------------------------------------------
% 步骤 1: 距离压缩 (Range Compression)
% -----------------------------------------------------------
fprintf('执行距离压缩...\n');

% 1.1 生成距离向时间轴
tau = (-Nrepl/2 : Nrepl/2-1) / Fr; % 脉冲持续时间内的采样点

% 1.2 生成理想的发射脉冲 (Reference Chirp)
% 文档中提到: Kr = 0.72135e12 Hz/s
ref_chirp = exp(1j * pi * Kr * tau.^2);

% 1.3 生成匹配滤波器 (频域)
% 为了提高效率，通常在频域做卷积。FFT 长度取 2 的幂次以加速
Nfft_r = 2^nextpow2(Nr + Nrepl); 
Ref_f = fft(ref_chirp, Nfft_r); 
Ref_f_conj = conj(Ref_f); % 匹配滤波器是信号共轭

% 1.4 对每一行进行压缩
rc_data = zeros(Na, Nr); % 初始化结果矩阵
for i = 1:Na
    % 对回波信号做 FFT
    echo_f = fft(raw_data_block(i, :), Nfft_r);
    
    % 频域相乘
    compressed_f = echo_f .* Ref_f_conj;
    
    % IFFT 回到时域
    compressed_t = ifft(compressed_f);
    
    % 截取有效数据 (消除卷积带来的边缘效应)
    % 这一步需要根据延迟仔细调整，通常取中间部分或开头部分
    % 这里做一个简单的对齐：
    rc_data(i, :) = compressed_t(1:Nr); 
end

% 显示距离压缩后的幅度图（看看有没有目标亮线）
figure;
imagesc(abs(rc_data));
title('距离压缩后 (Range Compressed)');
xlabel('距离向 (Range)'); ylabel('方位向 (Azimuth)');
colormap(gray); 
drawnow;

% -----------------------------------------------------------
% 步骤 2: 方位向 FFT (变换到 Range-Doppler 域)
% -----------------------------------------------------------
fprintf('变换到 R-D 域...\n');
rd_data = fft(rc_data, [], 1); % 对第一维(方位向)做 FFT

% -----------------------------------------------------------
% 步骤 3: 距离徙动校正 (RCMC)
% -----------------------------------------------------------
fprintf('执行 RCMC (简化版)...\n');
% 注意：RADARSAT 数据的 RCMC 比较复杂，需要计算精确的斜距 R0
% R(f_eta) = R0 * (1 + (lambda * f_eta / (2*Vr))^2 )
% 这里为了演示流程，暂时跳过插值 RCMC，直接做方位压缩
% 实际场景中，如果不做 RCMC，图像会有散焦，但目标能看出来。

% -----------------------------------------------------------
% 步骤 4: 方位压缩 (Azimuth Compression)
% -----------------------------------------------------------
fprintf('执行方位压缩...\n');

% 4.1 计算多普勒参数
lambda = c / f0;
Vr = 7062; % 等效雷达速度 (m/s)，需要根据卫星轨道参数估算，这里取典型值
% 方位频率轴
fa = linspace(-PRF/2, PRF/2, Na); 

% 4.2 计算最近距离 (Slant Range)
% 每个距离门的距离
range_axis = R0 + (0:Nr-1) * (c / (2 * Fr));

ac_data = zeros(Na, Nr); % 结果矩阵

% 对每个距离门(列)进行方位压缩
for j = 1:Nr
    R = range_axis(j);
    
    % 构造方位匹配滤波器 H_az(f_eta)
    % 相位 = 4 * pi * R * D(f_eta, Vr) / lambda
    % D(f_eta, Vr) = sqrt(1 - (lambda * fa / (2*Vr))^2)
    D = sqrt(1 - (lambda * fa / (2*Vr)).^2);
    H_az = exp(1j * 4 * pi * R * D / lambda);
    
    % 在 R-D 域相乘 (频域匹配滤波)
    % 注意：rd_data 的方位向 FFT 通常会有 fftshift 问题，需要对齐
    column_f = rd_data(:, j).';
    
    % 简单的匹配滤波：信号 * 滤波器的共轭
    compressed_az_f = column_f .* conj(H_az);
    
    % IFFT 回到图像域
    ac_data(:, j) = ifft(compressed_az_f);
end

% -----------------------------------------------------------
% 5. 显示最终结果
% -----------------------------------------------------------
figure;
imagesc(abs(ac_data));
title('最终成像结果 (Magnitude)');
xlabel('距离向 (Range)'); ylabel('方位向 (Azimuth)');
colormap(gray);
