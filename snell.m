function alive=snell(s_z,n_medium,k_medium,n_subs,k_subs)
%% ========================================================================
% Snell定律：判断光子是否透射或反射
% 功能：根据Fresnel反射率判断光子是否在界面反射或透射
% ========================================================================
% 输入参数：
%   s_z - 光子z方向分量（判断方向：>0向下，<0向上）
%   n_medium, k_medium - 介质（涂层）的折射率和消光系数
%   n_subs, k_subs - 基底的折射率和消光系数
% 输出参数：
%   alive - 1表示反射（继续在涂层内），0表示透射（离开涂层）
% ========================================================================

    cos_teta=abs(s_z);  % 入射角余弦
    
    % 判断光子方向，确定外部介质
    if (s_z>0)
       % 向下传播：外部介质为基底
       n_outside=n_subs;
       k_outside=k_subs; 
    else
       % 向上传播：外部介质为空气
       n_outside=1;
       k_outside=0;
    end
    
    % ========== 计算Fresnel反射率 ==========
    if (n_outside==n_medium && k_outside==k_medium)
        % 折射率相同：无反射，直接透射
        reflectance=0;
    elseif (cos_teta>0.9999)
        % 垂直入射：简化公式
        reflectance=((n_medium-n_outside)^2+(k_medium-k_outside)^2)/((n_medium+n_outside)^2+(k_medium+k_outside)^2);
    elseif (cos_teta<(0.0001))
        % 掠入射：全反射
        reflectance=1;
    else
        % 一般情况：使用Fresnel公式（考虑吸收介质）
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan2=(n_medium-1i*k_medium)/(n_outside-1i*k_outside);  % 相对折射率
        sin_x2=sin_teta*carpan2;
        cos_x2=sqrt(1-sin_x2*sin_x2);
        
        % 平行和垂直偏振分量的反射率
        E_parallel=(cos_teta/cos_x2-carpan2)/(cos_teta/cos_x2+carpan2);
        R_parallel=E_parallel*conj(E_parallel);
        E_orth=-(cos_x2/cos_teta-carpan2)/(cos_x2/cos_teta+carpan2);
        R_orth=E_orth*conj(E_orth);
        
        % 非偏振光：取平均
        reflectance=real(R_parallel+R_orth)*0.5;
    end
    
    % ========== 根据反射率判断是否反射 ==========
    Rastgele=rand();
    if Rastgele>reflectance
        % 透射：光子离开涂层
        alive=0;
    else
        % 反射：光子继续在涂层内传输
        alive=1;
    end
end
