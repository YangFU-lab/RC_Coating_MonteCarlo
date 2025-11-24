%% ========================================================================
% 辐射制冷涂层蒙特卡洛模拟程序
% 功能：利用蒙特卡洛算法求解辐射制冷涂层的光谱反射/透射/吸收率
%       针对多种材料、多种尺寸颗粒的复合涂层系统
% ========================================================================

clear; 
clc;
close all;

%% ========================================================================
% 第一部分：定义粒子1 - TiO2的尺寸分布
% ========================================================================

% TiO2粒子尺寸分布数据（半径nm, 比例%）
% 数据来源：R902+ TiO2粒径分布
size_TiO2 = [105.70947	0.36325
             122.41994	2.37436
             141.77199	7.2006
             164.1832	13.59286
             190.13713	18.41023
             220.19386	19.8254
             255.00194	17.63953
             295.31244	12.39709
             341.99518	6.19551
             396.05753	1.81029
             458.66595	0.19088];

r_TiO2 = size_TiO2(:,1) * 10^-9;  % 转换为米
V_TiO2 = 4*pi/3 * r_TiO2.^3;      % 粒子体积
p_TiO2 = size_TiO2(:,2);          % 数量比例
% 归一化体积分数分布
vp_TiO2 = V_TiO2.*p_TiO2 / (sum(V_TiO2.*p_TiO2,'all'));

% TiO2总体积分数（可设置多个值进行参数扫描）
fv_total_TiO2 = [0.2, 0.4, 0.6];

%% ========================================================================
% 第二部分：定义粒子2 - 空心玻璃球的尺寸分布
% ========================================================================

% 空心玻璃球尺寸分布数据（直径nm, 比例%）
size_Hglass = [5182	0.46
               5780.5	0.6
               6447.5	0.76
               7192	0.96
               8022.5	1.19
               8948.5	1.45
               9983.5	1.74
               11135	2.08
               12420	2.43
               13855	2.81
               15455	3.2
               17240	3.63
               19230	4.04
               21450	4.44
               23925	4.84
               26685	5.18
               29765	5.47
               33205	5.64
               37040	5.74
               41315	5.7
               46085	5.52
               51410	5.23
               57350	4.8
               63970	4.28
               71355	3.72
               79595	3.12
               88785	2.55
               98815	2.02
               110250	1.55
               123200	1.16
               137450	0.86
               153350	0.61
               171050	0.43
               190800	0.31
               212800	0.21
               237350	0.15];

r_Hglass = size_Hglass(:,1)/2 * 10^-9;      % 外半径（米）
ri_Hglass = r_Hglass - 2.5*10^-6;           % 内半径（米，壁厚2.5微米）
V_Hglass = 4*pi/3 * r_Hglass.^3;             % 外球体积
p_Hglass = size_Hglass(:,2);                % 数量比例
% 归一化体积分数分布
vp_Hglass = V_Hglass.*p_Hglass / (sum(V_Hglass.*p_Hglass,'all'));

% 空心玻璃球总体积分数（可设置多个值进行参数扫描）
fv_total_Hglass = [0.1, 0.2];

%% ========================================================================
% 第三部分：定义计算参数
% ========================================================================

% 蒙特卡洛模拟参数
photon_number = 10^2;  % 光子束数量（可增加以提高精度，但计算时间更长）

% 波长范围（单位：米）
lamda = (300:10:2500)' * 10^-9;  % 300-2500 nm（紫外到近红外）
% lamda = (2500:10:16000)' * 10^-9;  % 可选：中红外范围
lamda_um = lamda * 10^6;  % 转换为微米

% 入射角度
polar_angle = 0;  % 入射角（度），0度为垂直入射
% polar_angle = linspace(0, 89.99999, 30);  % 可选：多角度计算
polar_angle_rad = polar_angle * pi/180;

%% ========================================================================
% 第四部分：设置材料光学常数
% ========================================================================

% 介质（基体）光学常数 - PMMA
n_medium = PMMA_n(lamda);
k_medium = PMMA_k(lamda);

% TiO2粒子光学常数
n_TiO2 = TiO2_n(lamda);
k_TiO2 = TiO2_k(lamda);  % 使用纯TiO2的消光系数

% 空心玻璃球：外壳（SiO2）和内核（空气）
n_shell = sio2_n(lamda);  % 外壳折射率
k_shell = sio2_k(lamda);  % 外壳消光系数
n_core = ones(length(lamda), 1);   % 内核折射率（空气）
k_core = zeros(length(lamda), 1);  % 内核消光系数

% 基底光学常数（空气）
n_sub = ones(length(lamda), 1);
k_sub = zeros(length(lamda), 1);

%% ========================================================================
% 第五部分：涂层厚度设置
% ========================================================================

thickness = 200 * 10^-6;  % 涂层厚度（米），200微米

%% ========================================================================
% 第六部分：计算表面反射率
% ========================================================================

% 计算每个波长下的表面反射率和折射角
teta_prime = zeros(length(lamda), 1);
sur_reflection = zeros(length(lamda), 1);

for i = 1:length(lamda)
    % 计算折射角
    teta_prime(i) = F_fresnel_2(n_medium(i), k_medium(i), polar_angle_rad) * 180/pi;
    
    % 计算Fresnel反射率
    cos_teta = cosd(polar_angle);
    sin_teta = sqrt(1 - cos_teta*cos_teta);
    carpan2 = 1/(n_medium(i) - 1i*k_medium(i));
    sin_x2 = sin_teta * carpan2;
    cos_x2 = sqrt(1 - sin_x2*sin_x2);
    carpan1 = cos_teta / cos_x2;
    carpan3 = cos_x2 / cos_teta;
    
    % 平行和垂直偏振分量的反射率
    E_parallel = (carpan1 - carpan2) / (carpan1 + carpan2);
    R_parallel = E_parallel * conj(E_parallel);
    E_orth = (carpan3 - carpan2) / (carpan3 + carpan2);
    R_orth = E_orth * conj(E_orth);
    
    % 平均反射率（非偏振光）
    reflectance = real(R_parallel + R_orth) * 0.5;
    sur_reflection(i) = reflectance;
end

%% ========================================================================
% 第七部分：初始化结果矩阵
% ========================================================================

% 光谱反射率和吸收率
ref_lamda = zeros(length(lamda), length(fv_total_TiO2), length(fv_total_Hglass));
abs_lamda = zeros(length(lamda), length(fv_total_TiO2), length(fv_total_Hglass));

% 散射参数矩阵（用于蒙特卡洛计算）
mu_tot_arr_1 = zeros(length(lamda), length(r_TiO2) + length(r_Hglass));      % 总消光系数
scat_prob_arr_1 = zeros(length(lamda), length(r_TiO2) + length(r_Hglass));   % 散射概率
g_arr_1 = zeros(length(lamda), length(r_TiO2) + length(r_Hglass));          % 不对称参数

% 临时变量
r_no_1 = zeros(length(lamda), 1);
a_no_1 = zeros(length(lamda), 1);

% 太阳光谱
Solar = I_solar(lamda);

% 太阳反射率（加权平均）
Rtotal = zeros(length(fv_total_TiO2), length(fv_total_Hglass));

%% ========================================================================
% 第八部分：主计算循环 - 参数扫描
% ========================================================================

fprintf('开始蒙特卡洛模拟计算...\n');
fprintf('TiO2体积分数: %s\n', mat2str(fv_total_TiO2));
fprintf('空心玻璃球体积分数: %s\n', mat2str(fv_total_Hglass));
fprintf('波长数量: %d\n', length(lamda));
fprintf('光子束数量: %d\n\n', photon_number);

bar = waitbar(0, '计算进度...');

% 遍历所有TiO2体积分数
for m = 1:length(fv_total_TiO2)
    % 遍历所有空心玻璃球体积分数
    for n = 1:length(fv_total_Hglass)
        % 计算各尺寸粒子的体积分数
        fv_TiO2 = fv_total_TiO2(m) * vp_TiO2;
        fv_Hglass = fv_total_Hglass(n) * vp_Hglass;
        
        tic
        % 遍历所有波长
        for i = 1:length(lamda)
            % 更新进度条
            progress = ((m-1)*length(fv_total_Hglass)*length(lamda) + ...
                       (n-1)*length(lamda) + (i-1)) / ...
                       (length(fv_total_TiO2)*length(fv_total_Hglass)*length(lamda));
            waitbar(progress, bar, sprintf('计算进度: %.1f%%', progress*100));
            
            % ========== 计算TiO2粒子的散射参数 ==========
            for j = 1:length(r_TiO2)
                % Mie散射计算
                x_1 = 2*pi*n_medium(i)*r_TiO2(j)/lamda(i);
                m1 = (n_TiO2(i) + 1i*k_TiO2(i)) / n_medium(i);
                mie_result = Mie(m1, x_1);
                
                Qscat_1 = mie_result(2);  % 散射效率
                Qabs_1 = mie_result(3);   % 吸收效率
                Cscat_1 = pi*r_TiO2(j)^2 * Qscat_1;  % 散射截面
                Cabs_1 = pi*r_TiO2(j)^2 * Qabs_1;    % 吸收截面
                
                % 计算散射和吸收系数
                alfa_1 = fv_TiO2(j) * Cscat_1 / V_TiO2(j);  % 散射系数
                beta_1 = fv_TiO2(j) * Cabs_1 / V_TiO2(j);   % 吸收系数
                
                % 存储参数
                mu_tot_arr_1(i, j) = alfa_1 + beta_1;                    % 总消光系数
                scat_prob_arr_1(i, j) = alfa_1 / (alfa_1 + beta_1);     % 散射概率
                g_arr_1(i, j) = mie_result(5);                           % 不对称参数
            end
            
            % ========== 计算空心玻璃球的散射参数 ==========
            for j = 1:length(r_Hglass)
                % 核壳结构Mie散射计算
                x_2 = 2*pi*n_medium(i)*ri_Hglass(j)/lamda(i);  % 内球尺寸参数
                y_2 = 2*pi*n_medium(i)*r_Hglass(j)/lamda(i);   % 外球尺寸参数
                m1_2 = n_core(i) / n_medium(i);                 % 内核相对折射率
                m2_2 = (n_shell(i) + 1i*k_shell(i)) / n_medium(i);  % 外壳相对折射率
                
                miecoated_result = Miecoated(m1_2, m2_2, x_2, y_2, 1);
                
                Qscat_2 = miecoated_result(2);  % 散射效率
                Qabs_2 = miecoated_result(3);  % 吸收效率
                Cscat_2 = pi*r_Hglass(j)^2 * Qscat_2;  % 散射截面
                Cabs_2 = pi*r_Hglass(j)^2 * Qabs_2;    % 吸收截面
                
                % 计算散射和吸收系数
                alfa_2 = fv_Hglass(j) * Cscat_2 / V_Hglass(j);  % 散射系数
                beta_2 = fv_Hglass(j) * Cabs_2 / V_Hglass(j);   % 吸收系数
                
                % 存储参数（索引从TiO2之后开始）
                idx = j + length(r_TiO2);
                mu_tot_arr_1(i, idx) = alfa_2 + beta_2;                    % 总消光系数
                scat_prob_arr_1(i, idx) = alfa_2 / (alfa_2 + beta_2);     % 散射概率
                g_arr_1(i, idx) = miecoated_result(5);                    % 不对称参数
            end
            
            % 基体材料的吸收系数
            mu_matrix = 4*pi*k_medium(i) / lamda(i);
            
            % ========== 蒙特卡洛模拟 ==========
            % 调用蒙特卡洛函数计算反射率和吸收率
            [r_no_1(i), a_no_1(i)] = monte_carlo_single_layer_multipigment(...
                photon_number, ...                    % 光子束数量
                sur_reflection(i), ...                % 表面反射率
                cosd(teta_prime(i)), ...             % 折射角余弦
                thickness, ...                        % 涂层厚度
                scat_prob_arr_1(i, :), ...            % 各粒子散射概率
                mu_tot_arr_1(i, :), ...               % 各粒子总消光系数
                mu_matrix, ...                        % 基体吸收系数
                g_arr_1(i, :), ...                    % 各粒子不对称参数
                n_medium(i), k_medium(i), ...         % 介质光学常数
                n_sub(i), k_sub(i));                  % 基底光学常数
            
            % 存储结果
            ref_lamda(i, m, n) = r_no_1(i);  % 反射率
            abs_lamda(i, m, n) = a_no_1(i);  % 吸收率
        end
        
        elapsed_time = toc;
        fprintf('完成 fv_TiO2=%.2f, fv_Hglass=%.2f，耗时: %.2f秒\n', ...
                fv_total_TiO2(m), fv_total_Hglass(n), elapsed_time);
        
        % 计算太阳反射率（加权平均）
        Rtotal(m, n) = trapz(lamda, ref_lamda(:, m, n).*Solar) / trapz(lamda, Solar);
        
        % 绘制光谱反射率
        figure('Position', [100, 100, 800, 500]);
        plot(lamda*1e6, ref_lamda(:, m, n), 'LineWidth', 2);
        xlabel('波长 (\mum)', 'FontSize', 12);
        ylabel('反射率', 'FontSize', 12);
        title(sprintf('光谱反射率 (fv_{TiO2}=%.2f, fv_{Hglass}=%.2f)', ...
              fv_total_TiO2(m), fv_total_Hglass(n)), 'FontSize', 13);
        legend(sprintf('R_{solar}=%.3f', Rtotal(m, n)), 'Location', 'best', 'FontSize', 11);
        grid on;
        set(gca, 'FontSize', 11);
    end
end

delete(bar);

%% ========================================================================
% 第九部分：绘制参数扫描结果
% ========================================================================

% 绘制太阳反射率随体积分数的变化
[M, N] = meshgrid(fv_total_Hglass, fv_total_TiO2);
figure('Position', [100, 100, 800, 600]);
surf(M, N, Rtotal);
xlabel('空心玻璃球体积分数', 'FontSize', 12);
ylabel('TiO_2 体积分数', 'FontSize', 12);
zlabel('太阳反射率', 'FontSize', 12);
title('太阳反射率参数扫描', 'FontSize', 14, 'FontWeight', 'bold');
colorbar;
colormap('jet');
set(gca, 'FontSize', 11);

%% ========================================================================
% 程序结束
% 计算结果说明：
% - ref_lamda: 光谱反射率 [波长 × TiO2体积分数 × 空心玻璃球体积分数]
% - abs_lamda: 光谱吸收率 [波长 × TiO2体积分数 × 空心玻璃球体积分数]
% - Rtotal: 太阳反射率 [TiO2体积分数 × 空心玻璃球体积分数]
% ========================================================================

fprintf('\n计算完成！\n');

