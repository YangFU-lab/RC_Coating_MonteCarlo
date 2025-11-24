function result = Miecoated(m1,m2,x,y,opt)



% Mie Efficiencies of coated spheres for given complex refractive-index

% ratios m1=m1'+im1", m2= m2'+im2" of kernel and coating, resp.,

% and size parameters x=k0*a, y=k0*b where k0= wave number in ambient 

% medium, a,b= inner,outer sphere radius, using complex Mie Coefficients

% an and bn for n=1 to nmax,

% s. Bohren and Huffman (1983) BEWI:TDD122, p. 181-185,483.

% Result: m1',m1", m2',m2", x,y, efficiencies for extinction (qext), 

% scattering (qsca), absorption (qabs), backscattering (qb), 

% asymmetry parameter (asy=<costeta>) and (qratio=qb/qsca).

% opt selects the function "Miecoated_ab.." for an and bn, n=1 to nmax.

% Note that 0<=x<=y;   C. M�tzler, August 2002.



% ========== 特殊情况处理 ==========
if x==y
    % 内核半径等于外壳半径：退化为实心球
    result=Mie(m1,y);
elseif x==0
    % 内核半径为零：退化为实心球（外壳材料）
    result=Mie(m2,y);
elseif m1==m2
    % 内核和外壳折射率相同：退化为实心球
    result=Mie(m1,y);
elseif x>0
    % ========== 正常情况：核壳结构 ==========

    % 计算所需的最大模式数
    nmax=round(2+y+4*y.^(1/3));
    n1=nmax-1;
    n=(1:nmax);
    cn=2*n+1;  % 权重系数
    c1n=n.*(n+2)./(n+1);
    c2n=cn./n./(n+1);
    y2=y.*y;

    % 根据opt参数选择计算方法
    if opt==1
        f=Miecoated_ab1(m1,m2,x,y);
    elseif opt==2
        f=Miecoated_ab2(m1,m2,x,y);
    elseif opt==3
        f=Miecoated_ab3(m1,m2,x,y);
    end

    % 提取Mie系数的实部和虚部
    anp=(real(f(1,:)));  % an系数实部
    anpp=(imag(f(1,:))); % an系数虚部
    bnp=(real(f(2,:)));  % bn系数实部
    bnpp=(imag(f(2,:))); % bn系数虚部

    % 准备计算不对称参数所需的数组
    g1(1:4,nmax)=[0; 0; 0; 0];
    g1(1,1:n1)=anp(2:nmax);
    g1(2,1:n1)=anpp(2:nmax);
    g1(3,1:n1)=bnp(2:nmax);
    g1(4,1:n1)=bnpp(2:nmax);   

    % 计算消光效率
    dn=cn.*(anp+bnp);
    q=sum(dn);
    qext=2*q./y2;

    % 计算散射效率
    en=cn.*(anp.*anp+anpp.*anpp+bnp.*bnp+bnpp.*bnpp);
    q=sum(en);
    qsca=2*q./y2;

    % 计算吸收效率
    qabs=qext-qsca;

    % 计算后向散射效率
    fn=(f(1,:)-f(2,:)).*cn;
    gn=(-1).^n;
    f(3,:)=fn.*gn;
    q=sum(f(3,:));
    qb=q*q'./y2;

    % 计算不对称参数
    asy1=c1n.*(anp.*g1(1,:)+anpp.*g1(2,:)+bnp.*g1(3,:)+bnpp.*g1(4,:));
    asy2=c2n.*(anp.*bnp+anpp.*bnpp);
    asy=4/y2*sum(asy1+asy2)/qsca;

    % 计算后向散射比
    qratio=qb/qsca;

    % 返回结果
    result=[qext qsca qabs qb asy qratio];

end