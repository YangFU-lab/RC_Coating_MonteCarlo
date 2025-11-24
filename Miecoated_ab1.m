function result = Miecoated_ab1(m1,m2,x,y)
%% ========================================================================
% 核壳结构Mie系数计算（方法1）
% 功能：计算核壳结构球形粒子的Mie散射系数an和bn
% ========================================================================
% 输入参数：
%   m1 - 内核相对折射率（复数，m1'+im1"）
%   m2 - 外壳相对折射率（复数，m2'+im2"）
%   x - 内核尺寸参数（x = 2πn_medium*a/λ，a为内核半径）
%   y - 外壳尺寸参数（y = 2πn_medium*b/λ，b为外壳半径）
% 输出参数：
%   result - [an; bn]，第一行为an系数，第二行为bn系数
% ========================================================================
% 参考文献：Bohren & Huffman (1983) BEWI:TDD122, p. 483
%           使用递推关系(4.89)计算Dn函数
%           针对有损耗材料进行了数值优化，避免溢出和下溢
% ========================================================================

% ========== 计算参数 ==========
m=m2./m1;  % 外壳与内核的折射率比

% Bessel函数的参数
u=m1.*x;  % 内核参数
v=m2.*x;  % 外壳内表面参数
w=m2.*y;  % 外壳外表面参数

% 计算所需的最大模式数
nmax=round(2+y+4*y.^(1/3));
mx=max(abs(m1*y),abs(m2*y));
nmx=round(max(nmax,mx)+16);
nmax1=nmax-1;
n=(1:nmax);

% ========== 计算Dn函数（对三个参数u, v, w分别计算） ==========
% 根据Bohren & Huffman (1983) Eq. (4.89)

% Dn(u) - 内核参数
dnx(nmx)=0+0i;
z=u;
for j=nmx:-1:2
    dnx(j-1)=j./z-1/(dnx(j)+j./z);
end
dnu=dnx(n);

% Dn(v) - 外壳内表面参数
z=v;
for j=nmx:-1:2
    dnx(j-1)=j./z-1/(dnx(j)+j./z);
end
dnv=dnx(n);

% Dn(w) - 外壳外表面参数
z=w;
for j=nmx:-1:2
    dnx(j-1)=j./z-1/(dnx(j)+j./z);
end
dnw=dnx(n);

% ========== 计算Psi、Chi和Gsi函数及其导数 ==========
nu = (n+0.5);

% 计算各位置的球Bessel函数
sv= sqrt(0.5*pi*v);
pv= sv.*besselj(nu, v);  % 外壳内表面的Psi函数

sw= sqrt(0.5*pi*w);
pw= sw.*besselj(nu, w);  % 外壳外表面的Psi函数

sy= sqrt(0.5*pi*y);
py= sy.*besselj(nu, y);  % 环境介质中的Psi函数
p1y=[sin(y), py(1:nmax1)];  % Psi函数的导数

chv= -sv.*bessely(nu,v);  % 外壳内表面的Chi函数
chw= -sw.*bessely(nu,w);  % 外壳外表面的Chi函数
chy= -sy.*bessely(nu,y);  % 环境介质中的Chi函数
ch1y= [cos(y), chy(1:nmax1)];  % Chi函数的导数

gsy= py-i*chy;   % 环境介质中的Hankel函数
gs1y= p1y-i*ch1y;  % 其导数

% ========== 计算U、V、F函数（避免Riccati-Bessel函数乘积） ==========
uu=m.*dnu-dnv;  % U函数
vv=dnu./m-dnv;  % V函数

fv=pv./chv;  % F函数（外壳内表面）
fw=pw./chw;  % F函数（外壳外表面）

% 计算中间变量
ku1=uu.*fv./pw;
kv1=vv.*fv./pw;
ku2=uu.*(pw-chw.*fv)+(pw./pv)./chv;
kv2=vv.*(pw-chw.*fv)+(pw./pv)./chv;

dns1=ku1./ku2;  % Dn波浪号（第一部分）
gns1=kv1./kv2;  % Gn波浪号（第一部分）

% ========== 计算Dn_Schlange和Gn_Schlange ==========
dns=dns1+dnw;  % 完整的Dn波浪号
gns=gns1+dnw;  % 完整的Gn波浪号

% 计算最终系数所需的中间量
a1=dns./m2+n./y;
b1=m2.*gns+n./y;

% ========== 计算Mie系数 ==========
an=(py.*a1-p1y)./(gsy.*a1-gs1y);  % 电偶极模式系数
bn=(py.*b1-p1y)./(gsy.*b1-gs1y);  % 磁偶极模式系数

% 返回结果
result=[an; bn];
