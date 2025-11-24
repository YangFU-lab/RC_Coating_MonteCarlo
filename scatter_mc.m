function [s_x_,s_y_,s_z_]=scatter_mc(g,s_x,s_y,s_z)
%% ========================================================================
% 散射方向计算（Henyey-Greenstein相位函数）
% 功能：根据不对称参数g计算散射后的新方向向量
% ========================================================================
% 输入参数：
%   g - 不对称参数（-1到1，g>0前向散射，g<0后向散射，g=0各向同性）
%   s_x, s_y, s_z - 散射前的方向向量（归一化）
% 输出参数：
%   s_x_, s_y_, s_z_ - 散射后的方向向量（归一化）
% ========================================================================

    % ========== 步骤1：计算散射角（极角） ==========
    % 使用Henyey-Greenstein相位函数采样散射角
    rnd=rand();
    carpan=(1 - g*g)/(1 - g + 2*g*rnd);
    cos_theta=(1 + g*g - carpan*carpan)/(2*g);  % 散射角余弦
    sin_theta=sqrt(1 - cos_theta*cos_theta);     % 散射角正弦
    
    % ========== 步骤2：计算方位角 ==========
    R_azimuth=rand();
    phi=2*pi*R_azimuth;  % 方位角（0到2π，均匀分布）
    cos_phi=cos(phi);
    if phi<pi
        sin_phi=sqrt(1-cos_phi*cos_phi);
    else
        sin_phi=-sqrt(1-cos_phi*cos_phi);
    end
    
    % ========== 步骤3：计算新方向向量 ==========
    % 将散射角从局部坐标系转换到全局坐标系
    if (s_z==1)
        % 特殊情况：原方向沿+z轴
        s_x_ = sin_theta*cos_phi;
        s_y_ = sin_theta*sin_phi;
        s_z_ = cos_theta;                                           
    elseif (s_z==-1)
        % 特殊情况：原方向沿-z轴
        s_x_ = sin_theta*cos_phi;
        s_y_ = sin_theta*sin_phi;
        s_z_ = -cos_theta;
    else     
        % 一般情况：坐标旋转
        denom = sqrt(1 - s_z*s_z);
        s_x_ = sin_theta*(s_x * s_z * cos_phi - s_y * sin_phi) / denom + s_x * cos_theta;
        s_y_ = sin_theta*(s_y * s_z * cos_phi + s_x * sin_phi) / denom + s_y * cos_theta;
        s_z_ = -denom*sin_theta*cos_phi + s_z*cos_theta;
    end
end