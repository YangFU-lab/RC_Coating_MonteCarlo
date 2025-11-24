%% ========================================================================
% 多层梯度RC涂层光谱性能计算（单一颗粒材料）
% 功能：对于单一颗粒材料，计算具有浓度或尺寸梯度的辐射冷却涂层的光谱性能
%       采用分层计算思路，每层的体积分数和颗粒半径可以不同（梯度分布）
% ========================================================================
% 主要步骤：
%   1. 参数设置（波长、角度、厚度、层数等）
%   2. 定义梯度分布（体积分数和颗粒半径的线性变化）
%   3. 计算表面反射率
%   4. 对每个波长和每层计算Mie散射参数
%   5. 使用蒙特卡洛方法计算反射率和吸收率
%   6. 计算太阳反射率和红外吸收率
%   7. 绘制光谱结果
% ========================================================================

clc
clear all
close all

% ========== 第一部分：基本参数设置 ==========
% 蒙特卡洛模拟参数
photon_number = 10^5;  % 光子束数量

% 波长范围定义（单位：米）
lamda = [(300:10:2500), (2550:50:30000)]' * 10^-9;  % 全光谱范围
lamda_um = lamda * 10^6;  % 转换为微米
lamda_color = (360:830)' * 10^-9;  % 可见光范围（未使用）

% 入射角度设置
polar_angle = 0;  % 入射角（度），0度表示垂直入射
% polar_angle = linspace(0, 89.99999, 30);  % 可选：多角度计算
polar_angle_rad = polar_angle * pi / 180;  % 转换为弧度

% ========== 第二部分：光学常数定义 ==========
% 介质（基体）的光学常数
n_medium = PDMS101_n(lamda);  % 折射率
k_medium = PDMS101_k(lamda);  % 消光系数

% 颗粒（颜料）的光学常数
n_pigment = TiO2_n(lamda);  % 折射率
k_pigment = TiO2_k(lamda);  % 消光系数

% 基底的光学常数
n_sub = ones(length(lamda), 1);  % 折射率（空气）
k_sub = zeros(length(lamda), 1);  % 消光系数

% ========== 第三部分：涂层结构参数 ==========
% 厚度和层数
thickness = 100 * 10^-6;  % 总厚度（米）
thickness_um = thickness * 10^6;  % 转换为微米
layer_number = 5;  % 层数

% ========== 第四部分：梯度分布定义 ==========
% 体积分数梯度（从顶部到底部的线性变化）
f_v_top = 0.05;  % 顶部体积分数
f_v_bot = 0.05;  % 底部体积分数
f_v = zeros(layer_number, 1);  % 各层体积分数

% 颗粒半径梯度（从顶部到底部的线性变化）
pigment_r_top = 200 * 10^-9;  % 顶部颗粒半径（米）
pigment_r_bot = 200 * 10^-9;  % 底部颗粒半径（米）
pigment_r = zeros(layer_number, 1);  % 各层颗粒半径

% 计算各层的体积分数和颗粒半径（线性插值）
for ii = 1:layer_number
    f_v(ii) = f_v_top + (ii-1) * (f_v_bot - f_v_top) / layer_number;
    pigment_r(ii) = pigment_r_top + (ii-1) * (pigment_r_bot - pigment_r_top) / layer_number;
end

% 计算各层颗粒体积
V = 4 * pi * pigment_r.^3 / 3;  % 单个颗粒体积

% ========== 第五部分：表面反射率计算 ==========
% 计算入射界面的Fresnel反射率
teta_prime = zeros(length(lamda), 1);  % 折射角（度）
sur_reflection = zeros(length(lamda), 1);  % 表面反射率

for i = 1:length(lamda)
    % 计算折射角
    teta_prime(i) = F_fresnel_2(n_medium(i), k_medium(i), polar_angle_rad) * 180 / pi;
    
    % 计算Fresnel反射率（考虑平行和垂直偏振）
    cos_teta = cosd(polar_angle);
    sin_teta = sqrt(1 - cos_teta * cos_teta);
    carpan2 = 1 / (n_medium(i) - 1i * k_medium(i));
    sin_x2 = sin_teta * carpan2;
    cos_x2 = sqrt(1 - sin_x2 * sin_x2);
    carpan1 = cos_teta / cos_x2;
    carpan3 = cos_x2 / cos_teta;
    
    % 平行偏振反射率
    E_parallel = (carpan1 - carpan2) / (carpan1 + carpan2);
    R_parallel = E_parallel * conj(E_parallel);
    
    % 垂直偏振反射率
    E_orth = (carpan3 - carpan2) / (carpan3 + carpan2);
    R_orth = E_orth * conj(E_orth);
    
    % 平均反射率
    reflectance = real(R_parallel + R_orth) * 0.5;
    sur_reflection(i) = reflectance;
end

% ========== 第六部分：结果数组初始化 ==========
ref_lamda_1 = zeros(length(lamda), 1);  % 反射率光谱
abs_lamda_1 = zeros(length(lamda), 1);  % 吸收率光谱

% 各层的Mie散射参数数组
mu_tot_arr_1 = zeros(length(lamda), layer_number);  % 总消光系数
scat_prob_arr_1 = zeros(length(lamda), layer_number);  % 散射概率
g_arr_1 = zeros(length(lamda), layer_number);  % 不对称参数

% 蒙特卡洛结果
r_no_1 = zeros(length(lamda), 1);  % 反射光子数
a_no_1 = zeros(length(lamda), 1);  % 吸收光子数

% ========== 第七部分：太阳光谱和黑体辐射定义 ==========
Rtotal = zeros(1, 1);  % 太阳反射率
Atotal = zeros(1, 1);  % 红外吸收率

% 波长范围划分
lam_solar = lamda(1:221);  % 太阳光谱范围（0.3-2.5微米）
lam_IR = lamda(222:571);  % 红外光谱范围（2.55-30微米）

% 太阳光谱辐照度
Solar = I_solar(lam_solar);

% 黑体辐射（300K）
T = 300;  % 温度（K）
BB = I_bb(lam_IR, T);  % 黑体辐射光谱

% ========== 第八部分：主计算循环 ==========
% 对每个波长计算涂层的光谱性能
bar = waitbar(0, '计算中...');
tic

for i = 1:length(lamda)
    waitbar(i / length(lamda), bar);
    
    % 对每层计算Mie散射参数
    for j = 1:layer_number
        % 计算尺寸参数
        n_medium_sc = n_medium(i);
        x_1 = 2 * pi * n_medium_sc * pigment_r(j) / lamda(i);
        
        % 计算相对折射率
        m1 = (n_pigment(i) + 1i * k_pigment(i)) / n_medium_sc;
        
        % 计算Mie散射效率
        mie_result = Mie(m1, x_1);
        Qsca_1 = mie_result(2);  % 散射效率
        Qabs_1 = mie_result(3);  % 吸收效率
        
        % 计算散射和吸收截面
        Csca_1 = pi * pigment_r(j)^2 * Qsca_1;
        Cabs_1 = pi * pigment_r(j)^2 * Qabs_1;
        
        % 计算散射和吸收系数
        alfa_1 = f_v(j) * Csca_1 / V(j);  % 散射系数
        beta_1 = f_v(j) * Cabs_1 / V(j) + 4 * pi * k_medium(i) / lamda(i);  % 吸收系数（包含介质吸收）
        
        % 存储各层的参数
        mu_tot_arr_1(i, j) = alfa_1 + beta_1;  % 总消光系数
        scat_prob_arr_1(i, j) = alfa_1 / (alfa_1 + beta_1);  % 散射概率
        g_arr_1(i, j) = mie_result(5);  % 不对称参数
    end
    
    % 使用蒙特卡洛方法计算反射率和吸收率
    [r_no_1(i), a_no_1(i)] = monte_carlo_multi_layer_single_pigment(...
        photon_number, ...
        sur_reflection(i), ...
        cosd(teta_prime(i)), ...
        thickness, ...
        layer_number, ...
        scat_prob_arr_1(i, :), ...
        mu_tot_arr_1(i, :), ...
        g_arr_1(i, :), ...
        n_medium(i), ...
        k_medium(i), ...
        n_sub(i), ...
        k_sub(i));
    
    % 存储结果
    ref_lamda_1(i) = r_no_1(i);
    abs_lamda_1(i) = a_no_1(i);
end

close(bar);
toc

% ========== 第九部分：计算加权平均性能 ==========
% 计算太阳反射率（加权平均）
Rtotal = trapz(lam_solar, ref_lamda_1(1:221) .* Solar) / trapz(lam_solar, Solar);

% 计算红外吸收率（加权平均）
Atotal = trapz(lam_IR, abs_lamda_1(222:571) .* BB) / trapz(lam_IR, BB);

% ========== 第十部分：结果可视化 ==========
figure('Position', [100, 100, 1000, 600]);
plot(lamda * 1e6, ref_lamda_1, 'r-', 'LineWidth', 2);
hold on;
plot(lamda * 1e6, abs_lamda_1, 'b-', 'LineWidth', 2);
xlabel('波长 (\mum)', 'FontSize', 12);
ylabel('反射率/吸收率', 'FontSize', 12);
title('RC涂层光谱性能', 'FontSize', 14);
legend(['厚度=', num2str(thickness_um), ' \mum, 层数=', num2str(layer_number), ...
    ', R_{solar}=', num2str(Rtotal, '%.3f'), ', A_{IR}=', num2str(Atotal, '%.3f')], ...
    'FontSize', 11, 'Location', 'best');
grid on;
xlim([0.3, 30]);

% 显示计算结果
fprintf('\n计算结果：\n');
fprintf('太阳反射率 R_{solar} = %.4f\n', Rtotal);
fprintf('红外吸收率 A_{IR} = %.4f\n', Atotal);
fprintf('涂层厚度 = %.1f 微米\n', thickness_um);
fprintf('层数 = %d\n', layer_number);
