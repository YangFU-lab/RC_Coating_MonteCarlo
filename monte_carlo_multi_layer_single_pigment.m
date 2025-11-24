function [r_tot,a_tot]=monte_carlo_multi_layer_single_pigment(photon_number,s_ref,cos_gecen,h,layer_num,scat_prob_1,mu_tot_1,g_1,n_medium,k_medium,n_subs,k_subs)
%% ========================================================================
% 蒙特卡洛模拟：多层单颗粒系统
% 功能：使用蒙特卡洛方法模拟光子在多层涂层中的传输过程
%       每层包含单一类型的颗粒，不同层的颗粒参数可以不同
% ========================================================================
% 输入参数：
%   photon_number - 光子束数量（蒙特卡洛采样数）
%   s_ref - 表面反射率（入射界面的Fresnel反射率）
%   cos_gecen - 折射角余弦值（入射界面折射后的方向）
%   h - 涂层总厚度（米）
%   layer_num - 层数
%   scat_prob_1 - 各层的散射概率 [1×layer_num]
%   mu_tot_1 - 各层的总消光系数 [1×layer_num]
%   g_1 - 各层的不对称参数 [1×layer_num]
%   n_medium, k_medium - 介质（基体）的折射率和消光系数
%   n_subs, k_subs - 基底的折射率和消光系数
% 输出参数：
%   r_tot - 总反射率
%   a_tot - 总吸收率
%   t_tot - 总透射率（未返回但已计算）
% ========================================================================

% 初始化结果累加器
r_tot=0;  % 反射光子数
a_tot=0;  % 吸收光子数
t_tot=0;  % 透射光子数

% 对每个光子束进行蒙特卡洛模拟
for i=1:photon_number
    % 初始化单个光子的状态
    r_no=0;  % 是否反射
    a_no=0;  % 是否吸收
    t_no=0;  % 是否透射
    
    % ========== 步骤1：判断是否在表面反射 ==========
    Rastgele=rand();
    if Rastgele>s_ref
        % 光子进入涂层
        alive=1;
    else
        % 光子在表面反射
        alive=0;
        r_no=1;
    end
    
    % ========== 步骤2：初始化光子位置和方向 ==========
    x=0;
    y=0;
    z=0;  % 初始位置在涂层顶部
    s_x=0;
    s_y=sqrt(1-cos_gecen*cos_gecen);  % y方向分量
    s_z=cos_gecen;                    % z方向分量（垂直于涂层）
    
    % ========== 步骤3：设置多层参数 ==========
    mu_tot=mu_tot_1;  % 各层的总消光系数
    mu_tot_base=mu_tot_1(1);  % 基准层的消光系数
    mu_ratio=mu_tot/mu_tot_base;  % 各层相对于基准层的消光系数比
    % 计算自由程长度（基于基准层）
    l_beta=-log(rand())/mu_tot_base;
    
    % 每层厚度
    t_layer=h/layer_num;
    %layer_flag=ceil(z/t_layer);
    %layer_flag=1;
    
    % ========== 步骤4：光子传输循环 ==========
    while alive
        % 确定当前所在层
        if s_z>0
            % 向下传播
            layer_flag=fix(z/t_layer)+1;  % 当前层号
            % 计算到下一层界面的距离（考虑消光系数比）
            l_w=mu_ratio(layer_flag)*(t_layer*layer_flag-z)/s_z;
            
            if l_w<l_beta
                min_index=0;  % 先到达层界面
                min_l=l_w/mu_ratio(layer_flag);  % 实际距离
            else
                min_index=1;  % 先发生相互作用
                min_l=l_beta/mu_ratio(layer_flag);  % 实际距离
            end
        else
            % 向上传播
            layer_flag=ceil(z/t_layer);  % 当前层号
            % 计算到上一层界面的距离（考虑消光系数比）
            l_w=-mu_ratio(layer_flag)*(z-t_layer*(layer_flag-1))/s_z;
            
            if l_w<l_beta
                min_index=0;  % 先到达层界面
                min_l=l_w/mu_ratio(layer_flag);  % 实际距离
            else
                min_index=1;  % 先发生相互作用
                min_l=l_beta/mu_ratio(layer_flag);  % 实际距离
            end
        end
        
        % 避免距离为零导致的死循环
        if min_l==0
            min_l=l_beta*1e-4;
        end
        
        % 更新光子位置
        x=x+min_l*s_x;
        y=y+min_l*s_y;
        z=z+min_l*s_z;

        % ========== 情况1：到达层界面或边界 ==========
        if (min_index==0)
            if (abs(z-h/2)-h/2)==0
                % 到达涂层边界（顶部或底部）
                alive=snell(s_z,n_medium,k_medium,n_subs,k_subs);
                if alive==0
                    % 光子离开涂层
                    if s_z>0
                        t_no=1;  % 向下透射（进入基底）
                    else
                        r_no=1;  % 向上反射（返回空气）
                    end
                else
                    % 光子在边界反射，继续在涂层内传输
                    l_beta=l_beta-l_w;
                    s_z=-s_z;  % 反转方向
                end
            else
                % 到达层间界面（层内传输，不反射）
                l_beta=l_beta-l_w;  % 更新剩余自由程
            end
        % ========== 情况2：与颗粒发生相互作用 ==========
        else
            % 获取当前层的颗粒参数
            g_=g_1(layer_flag);  % 不对称参数
            scat_prob_=scat_prob_1(layer_flag);  % 散射概率
                   
            random_no=rand();
            if random_no<scat_prob_
                % 发生散射：更新光子方向
                [s_x,s_y,s_z]=scatter_mc(g_,s_x,s_y,s_z);
                % 重新计算自由程长度
                l_beta=-log(rand())/mu_tot_base;
            else
                % 被颗粒吸收
                alive=0;
                a_no=1;
            end
        end
    end
    
    % 累加该光子的结果
    r_tot=r_tot+r_no;
    a_tot=a_tot+a_no;
    t_tot=t_tot+t_no;
end

% 计算平均反射率、吸收率和透射率
r_tot=r_tot/photon_number;
a_tot=a_tot/photon_number;
t_tot=t_tot/photon_number;