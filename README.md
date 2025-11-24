[README.md](https://github.com/user-attachments/files/23719992/README.md)
# 辐射冷却涂层蒙特卡洛模拟程序

## 作者信息

- **姓名**: Yang FU
- **研究网站**: [https://yfu-research.com/](https://yfu-research.com/)

## 引用

如果您在研究中使用了本代码，请引用以下文献，十分感谢：

**1. Fu, Y., An, Y., Xu, Y., Dai, J. G., & Lei, D. (2022). Polymer coating with gradient‐dispersed dielectric nanoparticles for enhanced daytime radiative cooling. *EcoMat*, 4(2), e12169.**  
https://doi.org/10.1002/eom2.12169

**2. Yalçın, R. A., Blandre, E., Joulain, K., & Drévillon, J. (2020). Colored radiative cooling coatings with nanoparticles. *ACS Photonics*, 7(5), 1312-1322.**  
https://doi.org/10.1021/acsphotonics.0c00513


本仓库包含两个用于计算辐射冷却（Radiative Cooling, RC）涂层光谱性能的MATLAB程序，均采用蒙特卡洛（Monte Carlo）方法进行光子传输模拟。

## 程序概述

### 1. `RC_coating_MC_MultPigmentMultiSize.m`
**功能：** 单层多材料多尺寸颗粒涂层的光谱性能计算

**特点：**
- 支持多种颗粒材料（如TiO₂、空心玻璃球等）
- 每种材料支持多种尺寸分布
- 支持核壳结构颗粒（如空心玻璃球）
- 可进行参数扫描（不同体积分数组合）
- 适用于复合涂层系统的优化设计

**应用场景：** 研究不同材料组合和体积分数对涂层性能的影响

### 2. `RC_coating_MC_MultLayerMultiDistribution.m`
**功能：** 多层梯度涂层的光谱性能计算

**特点：**
- 单一颗粒材料
- 多层结构设计
- 支持体积分数和颗粒尺寸的梯度分布（线性变化）
- 每层参数可独立设置
- 适用于梯度涂层设计

**应用场景：** 研究梯度分布对涂层性能的影响，优化涂层结构

## 主要功能

- **光谱计算：** 计算涂层在0.3-30 μm范围内的反射率和吸收率
- **性能评估：** 计算太阳反射率（Rsolar）和红外吸收率（AIR）
- **参数优化：** 支持体积分数、颗粒尺寸等参数扫描
- **可视化：** 自动绘制光谱曲线和参数扫描结果

## 依赖函数

### Mie散射计算
- `Mie.m` - 球形粒子Mie散射效率计算
- `Mie_ab.m` - Mie系数计算
- `Miecoated.m` - 核壳结构Mie散射计算
- `Miecoated_ab1.m` - 核壳结构Mie系数计算

### 蒙特卡洛模拟
- `monte_carlo_single_layer_multipigment.m` - 单层多颗粒蒙特卡洛模拟
- `monte_carlo_multi_layer_single_pigment.m` - 多层单颗粒蒙特卡洛模拟
- `scatter_mc.m` - 散射方向计算
- `snell.m` - 折射计算

### 光学常数
- `PMMA_n.m`, `PMMA_k.m` - PMMA光学常数
- `PDMS101_n.m`, `PDMS101_k.m` - PDMS光学常数
- `TiO2_n.m`, `TiO2_k.m` - TiO₂光学常数
- `sio2_n.m`, `sio2_k.m` - SiO₂光学常数

### 辅助函数
- `F_fresnel_2.m` - Fresnel反射率计算
- `I_solar.m` - 太阳光谱辐照度
- `I_bb.m` - 黑体辐射光谱

## 使用方法

### 程序1：多材料多尺寸涂层
1. 修改颗粒材料的光学常数函数调用
2. 设置各材料的尺寸分布数据
3. 设置体积分数扫描范围
4. 运行程序，查看光谱结果和参数扫描图

### 程序2：多层梯度涂层
1. 设置涂层总厚度和层数
2. 定义顶部和底部的体积分数/颗粒半径
3. 程序自动计算各层的梯度分布
4. 运行程序，查看光谱结果和性能指标

## 输出结果

- **光谱反射率/吸收率：** 随波长的变化曲线
- **太阳反射率（Rsolar）：** 加权平均反射率
- **红外吸收率（AIR）：** 加权平均吸收率（程序2）
- **参数扫描图：** 性能随参数变化的3D图（程序1）

## 注意事项

- 蒙特卡洛模拟的计算精度与光子束数量成正比，但计算时间也会增加
- 建议根据计算资源调整 `photon_number` 参数
- 波长范围可根据需要修改，但需确保光学常数函数支持相应范围
- **请在使用本代码时引用相关文献**（见上方引用部分）

