# IPRMSOLVER

本仓库为数值最优化求解器预研究项目的代码整理版本，主要记录原始—对偶松弛内点法（Interior-Point Relaxation Method, IPRM）在二次规划问题（Quadratic Programming, QP）上的算法推导、程序实现与数值实验过程。

项目以开源二阶锥规划求解器 QOCO 的代码结构为参考，在QOCO源码基础上，将 IPRM 方法改写并实现为一个可用于求解 QP 问题的算法程序，同时使用标准测试集对算法性能和数值稳定性进行初步评估。

## 已完成工作

目前已完成的主要工作包括：

1. 研读 QOCO 求解器论文与源码，梳理其内点法求解框架，包括问题设置、初始化、KKT 系统构造、方向求解、迭代更新、残差计算和停止准则等核心模块。

2. 针对二次规划问题推导 IPRM-QP 的主要迭代公式，并基于 QOCO 的代码组织方式实现对应算法流程，包括初始化、KKT 系统构造、搜索方向求解、线搜索更新、残差计算与停止判断等步骤。

3. 编写 MATLAB 调用与测试脚本，使用 MAROS、NETLIB 等标准测试集中的算例对 IPRM-QP 进行批量测试，并与 QOCO 的相关结果进行初步比较。

4. 针对实验中出现的数值不稳定和求解精度不足问题，对算法公式进行了等价变形和数量级调整，并重写部分核心迭代步骤，以改善变量尺度、矩阵条件数和求解稳定性。

5. 在 QP 实现基础上，进一步尝试向更一般的二阶锥规划问题（SOCP）扩展，初步整理了 SOCP 形式下的问题设置、数据处理与测试脚本。

## 目录说明

当前仓库保留了项目推进过程中的若干阶段性目录，主要包括：

```text
.
├── qoco-test/
│   ├── qoco/                 # QOCO 源码及相关测试接口
│   ├── maros/                # MAROS 测试算例
│   ├── netlib/               # NETLIB 测试算例
│   ├── qoco.m                # MATLAB 调用接口
│   ├── run_qoco_lp.m          # QOCO 在线性规划测试集上的运行脚本
│   └── run_scale_maros_qoco.m # QOCO 在缩放后 MAROS 算例上的测试脚本
│
├── iprm-test-2-clear/
│   ├── iprm/                 # IPRM-QP 核心实现
│   ├── maros/                # MAROS 测试算例
│   ├── netlib/               # NETLIB 测试算例
│   ├── iprm.m                # MATLAB 调用接口
│   ├── iprm_mex.cpp          # MATLAB MEX 接口
│   ├── run_iprm_lp.m          # IPRM 在线性规划测试集上的运行脚本
│   ├── run_scale_maros.m      # 缩放后 MAROS 测试脚本
│   └── *.csv / *.png / *.txt # 实验结果与对比图
│
├── iprm-socp/
│   ├── iprm/                 # 面向 SOCP 扩展的 IPRM 实现尝试
│   ├── maros/
│   ├── netlib/
│   ├── ScaleData.m            # 数值缩放预处理函数
│   └── *_results*.csv / *.png # 不同正则化参数下的实验结果
│
├── iprm-socp-2/
│   ├── dimacs-sedumi/         # DIMACS / SeDuMi 相关测试数据
│   ├── iprm/
│   ├── maros/
│   ├── netlib/
│   ├── run_dimacs.m           # DIMACS 测试脚本
│   ├── run_scale_dimacs.m     # 缩放后 DIMACS 测试脚本
│   └── *_results*.csv / *.png # SOCP 扩展实验结果
│
├── iprm-test-0-2x2-sparse/    # 早期小规模稀疏测试版本
├── iprm-test-2/               # 早期 IPRM-QP 测试版本
└── README.md
```

其中，`qoco-test/` 保留了 QOCO 基线与测试接口，`iprm-test-2-clear/` 是清理后的 IPRM-QP 实现与测试目录，`iprm-socp/` 和 `iprm-socp-2/` 记录了后续向 SOCP 扩展及数值稳定性实验的研究。


## 说明

本项目代码主要用于个人科研项目展示和过程记录。MAROS、NETLIB、DIMACS / SeDuMi 等测试数据集的使用应遵循其原始来源和许可要求。
