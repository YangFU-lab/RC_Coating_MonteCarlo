function [r_tot,a_tot,t_tot]=monte_carlo_single_layer_multipigment(photon_number,s_ref,cos_gecen,h,scat_prob_1,mu_tot_1,mu_matrix,g_1,n_medium,k_medium,n_subs,k_subs)
%% ========================================================================
% 蒙特卡洛模拟：单层多颗粒系统
% 功能：使用蒙特卡洛方法模拟光子在含多种颗粒的涂层中的传输过程
%       计算反射率、吸收率和透射率
% ========================================================================
% 输入参数：
%   photon_number - 光子束数量（蒙特卡洛采样数，越大精度越高）
%   s_ref - 表面反射率（入射界面的Fresnel反射率）
%   cos_gecen - 折射角余弦值（入射界面折射后的方向）
%   h - 涂层厚度（米）
%   scat_prob_1 - 各颗粒的散射概率 [1×N]（N为颗粒类型总数）
%   mu_tot_1 - 各颗粒的总消光系数 [1×N]
%   mu_matrix - 基体材料的吸收系数
%   g_1 - 各颗粒的不对称参数 [1×N]（Henyey-Greenstein相位函数参数）
%   n_medium, k_medium - 介质（基体）的折射率和消光系数
%   n_subs, k_subs - 基底的折射率和消光系数
% 输出参数：
%   r_tot - 总反射率
%   a_tot - 总吸收率
%   t_tot - 总透射率
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
    % 光子初始位置（涂层顶部）
    x=0;
    y=0;
    z=0;
    % 光子方向向量（归一化）
    s_x=0;
    s_y=sqrt(1-cos_gecen*cos_gecen);  % y方向分量
    s_z=cos_gecen;                     % z方向分量（垂直于涂层）
    
    % 计算总消光系数和基体吸收概率
    mu_tot=sum(mu_tot_1)+mu_matrix;  % 总消光系数
    prob_matrix=mu_matrix/mu_tot;    % 基体吸收概率
    % 计算自由程长度（使用指数分布）
    l_beta=-log(rand())/mu_tot;      % 到下一次相互作用（散射或吸收）的距离
    % ========== 步骤3：光子传输循环 ==========
    while alive   
        % 计算到边界的距离
        if (s_z>0)
            l_w = (h - z)/s_z;  % 到下边界（基底）的距离
        else
            l_w = -z/s_z;       % 到上边界（空气）的距离
        end
        
        % 判断先到达边界还是先发生相互作用
        if l_w<l_beta
            min_index=1;  % 先到达边界
            min_l=l_w;
        else
            min_index=2;  % 先发生相互作用
            min_l=l_beta;
        end
        
        % 更新光子位置
        x=x+min_l*s_x;
        y=y+min_l*s_y;
        z=z+min_l*s_z;
        
        % ========== 情况1：光子到达边界 ==========
        if (min_index==1)
            % 使用Snell定律判断是否透射或反射
            alive=snell(s_z,n_medium,k_medium,n_subs,k_subs);
            if (alive==0)
                % 光子离开涂层
                if s_z>0
                    t_no=1;  % 向下透射（进入基底）
                else
                    r_no=1;  % 向上反射（返回空气）
                end
            else
                % 光子在边界反射，继续在涂层内传输
                l_beta=l_beta-l_w;  % 更新剩余自由程
                s_z=-s_z;           % 反转z方向
            end
        % ========== 情况2：光子发生相互作用 ==========
        else
            random_no=rand();
            if random_no<prob_matrix
                % 被基体材料吸收
                a_no=1;
                alive=0;
            else
                % 与颗粒发生相互作用
                % 根据各颗粒的消光系数比例，随机选择遇到的颗粒类型
                random_no_1=rand();
                temp=1;
                while (random_no_1>(sum(mu_tot_1(1:temp))/sum(mu_tot_1)))
                    temp=temp+1;
                end
                pigment_flag=temp;  % 确定遇到的颗粒类型
                
                % 获取该颗粒的散射参数
                g_=g_1(pigment_flag);           % 不对称参数
                scat_prob_=scat_prob_1(pigment_flag);  % 散射概率
                
                % 判断是散射还是吸收
                random_no_2=rand();
                if random_no_2<scat_prob_
                    % 发生散射：更新光子方向
                    [s_x,s_y,s_z]=scatter_mc(g_,s_x,s_y,s_z);
                    % 重新计算自由程长度
                    l_beta=-log(rand())/mu_tot;
                else
                    % 被颗粒吸收
                    alive=0;
                    a_no=1;
                end
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