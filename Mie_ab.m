function result = Mie_ab(m,x)
%% ========================================================================
% Mie系数计算
% 功能：计算球形粒子的Mie散射系数an和bn（n=1到nmax）
% ========================================================================
% 输入参数：
%   m - 相对折射率（复数，m'+im"），相对于周围介质
%   x - 尺寸参数（x = 2πn_medium*a/λ），a为粒子半径
% 输出参数：
%   result - [an; bn]，第一行为an系数，第二行为bn系数
% ========================================================================
% 参考文献：Bohren & Huffman (1983) BEWI:TDD122, Eq. (4.88)
%           使用递推关系(4.89)计算Dn函数
% ========================================================================

% ========== 计算参数 ==========
z=m.*x;  % 复数参数z = m*x

% 计算所需的最大模式数
nmax=round(2+x+4*x.^(1/3));  % Bohren and Huffman (1983) p477
nmx=round(max(nmax,abs(z))+16);  % 递推计算所需的更大模式数

n=(1:nmax);
nu = (n+0.5);  % 用于Bessel函数的阶数

% ========== 计算球Bessel函数 ==========
sx=sqrt(0.5*pi*x);
px=sx.*besselj(nu,x);  % Psi函数（球Bessel函数j_n）
p1x=[sin(x), px(1:nmax-1)];  % Psi函数的导数

chx=-sx.*bessely(nu,x);  % Chi函数（球Neumann函数y_n）
ch1x=[cos(x), chx(1:nmax-1)];  % Chi函数的导数

% 计算Hankel函数（用于外行波）
gsx=px-i*chx;   % 第一类Hankel函数
gs1x=p1x-i*ch1x;  % 其导数

% ========== 计算Dn函数（递推关系） ==========
dnx(nmx)=0+0i;
for j=nmx:-1:2
    % 根据Bohren & Huffman (1983) Eq. (4.89)计算Dn(z)
    dnx(j-1)=j./z-1/(dnx(j)+j./z);
end
dn=dnx(n);  % Dn(z), n=1 to nmax

% ========== 计算中间变量 ==========
da=dn./m+n./x;  % 用于计算an的中间量
db=m.*dn+n./x;  % 用于计算bn的中间量

% ========== 计算Mie系数 ==========
% 根据Bohren & Huffman (1983) Eq. (4.88)
an=(da.*px-p1x)./(da.*gsx-gs1x);  % 电偶极模式系数
bn=(db.*px-p1x)./(db.*gsx-gs1x);  % 磁偶极模式系数

% 返回结果
result=[an; bn];
