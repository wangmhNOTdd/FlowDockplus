# PhysFlowDock: 基于双向协同投影的物理一致性全柔性蛋白质-配体对接

## 完整 Idea 描述（聚焦 Step 1 & Step 2，工程可行版本）

---

## 一、核心定位：我们与 CCFM 的本质区别

### 1.1 CCFM 做了什么？

经过仔细阅读 CCFM 原文（附录 B.1），其在分子对接任务中定义了四类约束：

- **C1 — 配体键长 (Bond Lengths)**
- **C2 — 配体键角 (Bond Angles)**
- **C3 — 配体内部碰撞排除 (Internal Steric Clash)**
- **C4 — 蛋白-配体最小距离/碰撞排除 (Minimum Protein–Ligand Contact)**

CCFM 在 FlexDock 的推理过程中嵌入了上述四类约束的投影算子。**CCFM 确实已经处理了蛋白-配体之间的碰撞问题**（Figure 2 是铁证：CCFM 将蛋白碳原子 CG 与配体氧原子 O1 之间的距离从 2.4Å 推至 3.8Å）。

### 1.2 CCFM 的关键假设与局限

根据 CCFM 附录 B.7.1 的实现细节，其投影算子的核心假设是：

> **蛋白质是静态刚性约束源（Rigid Constraint Source）。**

具体表现为：

1. **投影仅更新配体坐标** $x^L$：当发生蛋白-配体碰撞时，CCFM 只计算残差对配体坐标的雅可比 $\partial r / \partial x^L$，投影修正量 $\Delta x$ 仅施加于配体。蛋白原子 $x^P$ 在投影过程中**纹丝不动**。
2. **不检查蛋白内部合法性**：CCFM 对蛋白的键长、键角、内部碰撞没有任何约束。如果生成模型预测的蛋白构象本身存在物理违规，CCFM 无能为力。
3. **骨干模型限制**：CCFM 基于 FlexDock，后者使用"流形对接流（Manifold Docking Flow）"参数化配体刚体运动和蛋白侧链扭转角，本质上是**低维参数空间**的生成。蛋白的柔性仅限于侧链扭转角，而非全原子笛卡尔坐标。

### 1.3 PhysFlowDock 的核心创新

PhysFlowDock 基于 FlowDock——一个在**全原子笛卡尔坐标空间**同时生成蛋白和配体构象的全柔性模型。这意味着：

- FlowDock 的蛋白在生成过程中每个原子都在运动（不仅仅是侧链扭转）
- 蛋白的运动幅度可能很大（从 ESMFold 的 apo 态到 holo 态的全局形变）
- 蛋白和配体的碰撞问题需要**双向解决**

**PhysFlowDock 的核心术语：双向协同投影（Bilateral Cooperative Projection, BCP）**

| 特性 | CCFM + FlexDock | **PhysFlowDock + FlowDock** |
|:---|:---|:---|
| 蛋白表示 | 低维参数空间（扭转角） | **全原子笛卡尔坐标** |
| 碰撞投影对象 | 仅更新 $x^L$ | **同时更新 $x^L$ 和 $x^P$** |
| 蛋白在投影中的角色 | 刚性墙（不可移动） | **柔性参与者（可协同避让）** |
| 物理含义 | 配体被"推开" | **诱导契合（Induced-fit）：蛋白残基也为配体"让位"** |
| ODE 求解器 | 标准 Euler | **VD-ODE（需专门适配 CCFM 理论）** |

---

## 二、前置知识

### 2.1 FlowDock 的 VD-ODE 求解器

FlowDock 学习一个条件向量场 $v_\theta(x_t, t)$，但**不使用标准 Euler 积分**，而是使用 Variance Diminishing ODE (VD-ODE)：

$$
x_{n+1} = \underbrace{\text{clamp}\left(\frac{1-s}{1-t}\right)^\eta}_{\alpha_n} \cdot x_n + \underbrace{\text{clamp}\left(1 - \left(\frac{1-s}{1-t}\right)^\eta\right)}_{1 - \alpha_n} \cdot v_\theta(x_n, t_n)
$$

其中：
- $t_n = n/N$，$s = t_{n+1} = (n+1)/N$，$N$ 为总步数（默认 40）
- $\eta$ 为方差衰减参数
- $\text{clamp}$ 将值限制在 $[0, 1]$ 内
- $x_n = [x_n^P, x_n^L]$ 为蛋白和配体的联合坐标

**VD-ODE 的物理直觉**：在 $t$ 较小时（早期），$\alpha_n$ 接近 1，更新主要保留当前状态；在 $t$ 接近 1 时（晚期），$\alpha_n$ 接近 0，更新收敛到网络预测值 $v_\theta$。这意味着**早期运动幅度大，晚期精细调整**。

### 2.2 CCFM 的核心理论

CCFM 基于 OT 路径的线性关系 $x_t = (1-t)x_0 + tx_1$，推导出：在中间时刻 $t$ 对含噪样本 $x_t$ 的约束可以松弛为概率约束：

$$
\Pr_\xi\left[g(t^{-1}x_t - \xi) \leq 0\right] \geq 1 - \alpha
$$

其中 $\xi = (1-t)t^{-1}x_0$。对于高斯噪声源 $x_0 \sim \mathcal{N}(0, \sigma_0^2 I)$，该概率约束可以转化为**确定性的时间依赖容差**：

$$
g(x) \leq \beta(t)
$$

其中容差边界 $\beta(t)$ 随时间单调递减，在 $t=1$ 时降为 0（约束变为硬约束）。

**投影算子 $P_{cc}$**：

$$
x' = x - J^\top(JJ^\top + \lambda I)^{-1} \text{ReLU}(r - \beta(t))
$$

其中 $r$ 为约束残差向量，$J$ 为残差对坐标的雅可比矩阵。

---

## 三、PhysFlowDock 方法详述

### 3.1 总体架构：CCFM 投影嵌入 VD-ODE

**核心修改点：FlowDock 推理脚本中的 ODE 采样循环。**

原始 FlowDock 的采样循环（伪代码）：

```python
# FlowDock 原始推理 (inference.py 中的 sample 函数)
x = x_0  # 初始状态 [x_0^P, x_0^L]
for n in range(N):
    t = n / N
    s = (n + 1) / N
    v = model.forward(x, t)  # 神经网络前向传播
    alpha = clamp((1 - s) / (1 - t)) ** eta
    x = alpha * x + (1 - alpha) * v  # VD-ODE 更新
# x 即为最终预测的复合物结构
```

**PhysFlowDock 修改后的采样循环**：

```python
# PhysFlowDock 推理 (修改后)
x = x_0  # 初始状态 [x_0^P, x_0^L]
for n in range(N):
    t = n / N
    s = (n + 1) / N
    v = model.forward(x, t)  # 神经网络前向传播（不修改）
    alpha = clamp((1 - s) / (1 - t)) ** eta
    x_tilde = alpha * x + (1 - alpha) * v  # VD-ODE 更新（不修改）

    # ===== 以下为 PhysFlowDock 新增 =====
    # Step 1: 配体内部约束投影 (C1-C3) —— 仅更新 x^L
    x_tilde = ligand_internal_projection(x_tilde, phi(s), s)

    # Step 2: 蛋白-配体双向协同投影 (C4-C5) —— 同时更新 x^L 和 x^P
    x_tilde = bilateral_clash_projection(x_tilde, phi(s), s)

    x = x_tilde
# x 即为物理一致的复合物结构，无需 OpenMM 后处理
```

**关键设计决策：投影使用 $s = t_{n+1}$ 而非 $t_n$。** 因为 $x_\text{tilde}$ 是 $t_{n+1}$ 时刻的预测值，其容差应当按 $t_{n+1}$ 计算。

### 3.2 VD-ODE 与 CCFM 理论的兼容性分析

**问题**：CCFM 的理论推导基于标准 OT 路径 $x_t = (1-t)x_0 + tx_1$，对应标准 Euler 积分。FlowDock 使用的 VD-ODE 是一种非标准求解器，CCFM 理论能否直接适用？

**分析**：

1. **VD-ODE 的本质**：VD-ODE 可以理解为一种"加权 Euler"更新。当 $\eta = 1$ 时，VD-ODE 退化为标准线性插值（与 OT 路径一致）。当 $\eta > 1$ 时，路径在早期保留更多噪声，晚期更快收敛。
2. **CCFM 投影的适用性**：CCFM 投影算子 $P_{cc}$ 的本质是一个**后验修正步（Post-hoc Correction）**。无论 ODE 求解器使用何种格式，在每一步完成后，我们都可以对输出进行投影。关键是容差 $\beta(t)$ 的设计。
3. **容差调度的适配**：在 VD-ODE 下，$x_t$ 中的噪声成分不严格等于 $(1-t)\sigma_0$，而是受 $\eta$ 调制。因此，我们**不使用 CCFM 原文的解析容差公式**，而是采用一个**经验概率调度器** $\phi(t)$ 来控制约束松紧。

**结论**：CCFM 投影作为"即插即用"的后处理步骤，可以安全地嵌入 VD-ODE 循环。理论上的最优性可能因非标准路径而略有损失，但实践中完全可行。这也是一个重要的工程贡献——**首次将 CCFM 从标准 Euler FM 推广到 VD-ODE 求解器**。

### 3.3 概率调度器 $\phi(t)$ 的设计

我们不直接使用 CCFM 原文的 $\phi(t) = (t/2)^n$（该公式依赖于 OT 路径的精确高斯噪声结构），而是设计一个适配 VD-ODE 的经验调度器：

$$
\phi(t) = \begin{cases}
\phi_{\min} + (\phi_{\text{mid}} - \phi_{\min}) \cdot (t / t_s)^{n_1}, & t \leq t_s \\[6pt]
\phi_{\text{mid}} + (1 - \phi_{\text{mid}}) \cdot \left(\frac{t - t_s}{1 - t_s}\right)^{n_2}, & t > t_s
\end{cases}
$$

其中：
- $\phi_{\min} = 0.0$：$t=0$ 时不施加任何约束（纯噪声状态，约束无意义）
- $\phi_{\text{mid}} = 0.5$：切换点处的概率水平
- $t_s = 0.5$：切换时间点
- $n_1 = 0.5$（早期阶段，约束缓慢收紧）
- $n_2 = 0.1$（晚期阶段，约束迅速收紧至 $\phi \to 1$）

**当 $\phi(t)$ 确定后，容差边界为**：

$$
\beta_k(t) = (1 - \phi(t)) \cdot \beta_k^{\max}
$$

其中 $\beta_k^{\max}$ 为第 $k$ 类约束的最大容差（如键长约束 $\beta^{\max} = 0.5\text{Å}$，碰撞约束 $\beta^{\max} = 1.0\text{Å}$）。

当 $t \to 1$，$\phi(t) \to 1$，$\beta_k \to 0$，约束变为硬约束。

**设计理由**：
- **早期（$t < 0.5$）**：配体正在从随机位置向口袋大范围移动，过严的碰撞约束会阻碍配体穿越蛋白表面"探索"结合路径。
- **晚期（$t > 0.5$）**：分子已经接近最终结合态，必须严格禁止任何物理违规。

### 3.4 约束体系定义

#### 3.4.1 第一层约束：配体内部几何约束（Step 1）

**数据预处理**：从配体 SMILES 中使用 RDKit 解析分子图，提取所有共价键对、键角三元组、原子范德华半径。

**C1 — 共价键长约束**

对于每对键合原子 $(i, j) \in \mathcal{E}_{\text{bond}}$：

$$
g_{\text{bond,upper}}^{(i,j)}(x) = \|x_i^L - x_j^L\| - (1+\delta_b) \cdot \ell_{ij}^{\text{up}} \leq 0
$$
$$
g_{\text{bond,lower}}^{(i,j)}(x) = (1-\delta_b) \cdot \ell_{ij}^{\text{low}} - \|x_i^L - x_j^L\| \leq 0
$$

其中 $\delta_b = 0.25$（与 CCFM 一致），$\ell_{ij}^{\text{low/up}}$ 来自 RDKit 化学知识库。

**C2 — 键角约束**

对于每个键角三元组 $(i, j, k) \in \mathcal{E}_\angle$：

$$
g_{\angle,\text{upper}}^{(i,j,k)}(x) = \angle(x_i^L, x_j^L, x_k^L) - (1+\delta_\alpha) \cdot \alpha_{ijk}^{\text{up}} \leq 0
$$
$$
g_{\angle,\text{lower}}^{(i,j,k)}(x) = (1-\delta_\alpha) \cdot \alpha_{ijk}^{\text{low}} - \angle(x_i^L, x_j^L, x_k^L) \leq 0
$$

其中 $\delta_\alpha = 0.25$。

**C3 — 配体内部空间碰撞排除**

对于所有**非键合、非成角原子对**（拓扑距离 $\geq 3$）：

$$
g_{\text{intra}}^{(i,j)}(x) = (\rho_i + \rho_j) - \delta_s - \|x_i^L - x_j^L\| \leq 0
$$

其中 $\delta_s = 0.75\text{Å}$（与 CCFM 一致）。

#### 3.4.2 第二层约束：蛋白-配体界面约束（Step 2）

**C4 — 蛋白-配体范德华碰撞排除（核心约束）**

对于所有蛋白-配体原子对 $(i \in L, j \in P)$：

$$
g_{\text{clash}}^{(i,j)}(x) = (\rho_i + \rho_j) - \delta_{pl} - \|x_i^L - x_j^P\| \leq 0
$$

其中 $\delta_{pl} = 0.75\text{Å}$（与 CCFM 一致）。

**C5 — 结合口袋接触约束（仅在 $t > 0.7$ 时启用，可选）**

$$
g_{\text{contact}}(x) = d_{\min}(X^L, X^P) - d_{\text{contact}} \leq 0
$$

其中 $d_{\text{contact}} = 5.0\text{Å}$，防止配体在投影中被推离口袋过远。

---

### 3.5 投影算子的具体实现

#### 3.5.1 第一层投影：配体内部约束（仅更新 $x^L$）

此层与 CCFM 完全一致，**仅修改配体原子坐标**。

```python
def ligand_internal_projection(x, phi_t, t, max_iter=3, lam=1e-6):
    """
    第一层投影：修正配体内部键长、键角、自碰撞。
    仅更新 x 中的配体部分。

    Args:
        x: (N_total, 3) 全部原子坐标 [蛋白; 配体]
        phi_t: 当前概率水平 φ(t)
        t: 当前采样时间
        max_iter: Gauss-Newton 迭代次数
        lam: Tikhonov 正则化
    Returns:
        x: 投影后的坐标
    """
    x_P = x[:N_prot]        # 蛋白坐标（固定不动）
    x_L = x[N_prot:].clone() # 配体坐标（待修正）

    for _ in range(max_iter):
        x_L.requires_grad_(True)

        # 计算 C1, C2, C3 的残差
        r = compute_ligand_residuals(x_L, bond_pairs, bond_bounds,
                                      angle_triples, angle_bounds,
                                      nonbond_pairs, vdw_radii)

        # 应用时间依赖容差
        beta = (1 - phi_t) * beta_max_ligand  # 标量或向量
        r_active = torch.relu(r - beta)

        if r_active.max() < 1e-4:
            break  # 所有约束已满足

        # 筛选活跃约束
        active_mask = r_active > 0
        r_act = r_active[active_mask]

        # 计算雅可比矩阵 (仅对 x_L)
        J = torch.autograd.functional.jacobian(
            lambda z: compute_ligand_residuals(z, ...)[active_mask],
            x_L
        )  # shape: (|A|, N_lig * 3)

        J = J.reshape(r_act.shape[0], -1)

        # Gauss-Newton 更新
        delta = -J.T @ torch.linalg.solve(
            J @ J.T + lam * torch.eye(J.shape[0], device=J.device),
            r_act
        )
        x_L = (x_L.detach().reshape(-1) + delta).reshape(-1, 3)

    x = torch.cat([x_P, x_L.detach()], dim=0)
    return x
```

#### 3.5.2 第二层投影：双向协同投影（BCP，同时更新 $x^L$ 和 $x^P$）

**这是 PhysFlowDock 的核心创新。**

CCFM 在计算蛋白-配体碰撞投影时，仅计算：
$$
\Delta x^L = -J_L^\top(J_L J_L^\top + \lambda I)^{-1} r
$$

PhysFlowDock 计算完整的双向雅可比：
$$
J = \left[\frac{\partial r}{\partial x^L}, \frac{\partial r}{\partial x^P}\right]
$$

然后联合更新：
$$
\Delta \begin{bmatrix} x^L \\ x^P \end{bmatrix} = -J^\top(JJ^\top + \lambda I)^{-1} r
$$

**但直接这样做会带来一个问题**：蛋白有数千个原子，雅可比矩阵巨大。

**解决方案：局部活跃集 + 加权投影。**

我们只允许蛋白中**与碰撞相关的局部残基**参与投影（而非全部蛋白原子）：

```python
def bilateral_clash_projection(x, phi_t, t, max_iter=2, lam=1e-6, w_P=0.3):
    """
    第二层投影：蛋白-配体碰撞的双向协同投影。
    同时更新配体和蛋白的局部原子坐标。

    Args:
        x: (N_total, 3) 全部原子坐标
        phi_t: 当前概率水平
        t: 当前采样时间
        max_iter: Gauss-Newton 迭代次数
        lam: Tikhonov 正则化
        w_P: 蛋白位移权重 (0 < w_P ≤ 1)。
             w_P = 0 退化为 CCFM（仅动配体）。
             w_P = 1 表示蛋白和配体等权移动。
             默认 w_P = 0.3 表示蛋白承担 30% 的避让位移。
    Returns:
        x: 投影后的坐标
    """
    x_P = x[:N_prot].clone()  # 蛋白坐标
    x_L = x[N_prot:].clone()  # 配体坐标

    for _ in range(max_iter):
        # Step 1: 快速筛选碰撞候选对
        dist = torch.cdist(x_L, x_P)                        # (N_L, N_P)
        vdw_sum = vdw_L.unsqueeze(1) + vdw_P.unsqueeze(0)   # (N_L, N_P)
        threshold = vdw_sum - delta_pl                        # 碰撞阈值
        cutoff = threshold + 0.5  # 额外 0.5Å 缓冲区
        candidate_mask = dist < cutoff                        # 候选碰撞对

        if not candidate_mask.any():
            break  # 无碰撞

        # Step 2: 计算残差
        r_raw = threshold[candidate_mask] - dist[candidate_mask]  # >0 表示碰撞
        beta = (1 - phi_t) * beta_max_clash
        r_active = torch.relu(r_raw - beta)

        if r_active.max() < 1e-4:
            break

        active_in_candidates = r_active > 0
        r_act = r_active[active_in_candidates]

        # Step 3: 确定参与投影的蛋白原子集
        # 只有与碰撞直接相关的蛋白原子参与更新
        lig_indices, prot_indices = candidate_mask.nonzero(as_tuple=True)
        active_lig = lig_indices[active_in_candidates]
        active_prot = prot_indices[active_in_candidates]
        unique_prot_atoms = active_prot.unique()  # 涉及碰撞的蛋白原子

        # Step 4: 构建联合坐标向量 z = [x_L, x_P_active]
        x_P_active = x_P[unique_prot_atoms]  # 仅参与碰撞的蛋白原子
        z = torch.cat([x_L.reshape(-1), x_P_active.reshape(-1)])
        z.requires_grad_(True)

        # Step 5: 计算雅可比（对联合坐标 z）
        def compute_clash_residuals_joint(z_flat):
            z_L = z_flat[:N_lig*3].reshape(-1, 3)
            z_P_active = z_flat[N_lig*3:].reshape(-1, 3)
            # 构建完整蛋白坐标（大部分不变，仅活跃原子更新）
            z_P_full = x_P.clone()
            z_P_full[unique_prot_atoms] = z_P_active
            d = torch.cdist(z_L, z_P_full)
            r = threshold[candidate_mask] - d[candidate_mask]
            return r[active_in_candidates]

        J = torch.autograd.functional.jacobian(
            compute_clash_residuals_joint, z
        )  # shape: (|A|, dim_z)

        J = J.reshape(r_act.shape[0], -1)

        # Step 6: 加权 Gauss-Newton 更新
        # 构建权重矩阵 W：配体部分权重 1，蛋白部分权重 w_P
        dim_L = N_lig * 3
        dim_P_active = unique_prot_atoms.shape[0] * 3
        W = torch.ones(dim_L + dim_P_active, device=z.device)
        W[dim_L:] = w_P  # 蛋白部分的权重

        # 加权雅可比
        J_w = J * W.unsqueeze(0)  # 等效于 J @ diag(W)
        delta = -W * (J_w.T @ torch.linalg.solve(
            J_w @ J_w.T + lam * torch.eye(J.shape[0], device=J.device),
            r_act
        ))

        # Step 7: 应用更新
        z_new = z.detach() + delta
        x_L = z_new[:dim_L].reshape(-1, 3)
        x_P[unique_prot_atoms] = z_new[dim_L:].reshape(-1, 3)

    x = torch.cat([x_P.detach(), x_L.detach()], dim=0)
    return x
```

**权重参数 $w_P$ 的物理意义**：

- $w_P = 0$：退化为 CCFM 方案（蛋白不动，仅配体避让）
- $w_P = 0.3$（推荐默认值）：蛋白承担 30% 的避让位移。这模拟了真实的诱导契合效应——配体进入口袋时，蛋白侧链会轻微调整以容纳配体。
- $w_P = 1.0$：蛋白和配体等权移动。
- $w_P$ 可以设计为时间依赖的：$w_P(t) = w_P^{\max} \cdot t$（早期蛋白少动，晚期蛋白参与更多）

### 3.6 雅可比矩阵的高效计算策略

**方案选择**：使用 `torch.autograd` 自动微分。

**理由**：
1. 碰撞约束的雅可比矩阵是稀疏的（每个约束只涉及一对原子的 6 个坐标分量），手动解析虽然更快，但开发时间长且容易出错。
2. `torch.autograd.functional.jacobian` 在活跃集较小时（通常 $|\mathcal{A}| < 50$）计算速度完全可接受。
3. 在工程验证阶段使用自动微分，后续可以替换为手动解析以加速。

**手动解析雅可比（备选，用于后期优化）**：

对于碰撞约束 $g = (\rho_i + \rho_j) - \delta - \|x_i - x_j\|$：

$$
\frac{\partial g}{\partial x_i} = -\frac{x_i - x_j}{\|x_i - x_j\|}, \quad \frac{\partial g}{\partial x_j} = \frac{x_i - x_j}{\|x_i - x_j\|}
$$

每个约束的雅可比行是一个极稀疏向量（仅 6 个非零元素），可以直接构造稀疏矩阵。

### 3.7 计算复杂度分析

| 操作 | 数据规模 | 复杂度 | 估计耗时/步 |
|:---|:---|:---|:---|
| 配对距离矩阵 `torch.cdist` | $N_L \times N_P$（典型 40×2000） | GPU 并行 $O(N_L N_P)$ | ~1ms |
| 第一层 Jacobian + 投影 | 活跃约束数 $M_1$（典型 5-20） | $O(M_1^2 \cdot 3N_L)$ | ~1ms |
| 第二层 Jacobian + 投影 | 活跃约束数 $M_2$（典型 0-30） | $O(M_2^2 \cdot 3(N_L + N_P^{\text{active}}))$ | ~3ms |
| **每步总投影开销** | — | — | **~5ms** |
| **40 步总开销** | — | — | **~0.2s** |

FlowDock 的原始推理时间为 ~39s（神经网络前向 × 40 步），PhysFlowDock 的额外开销 ~0.2s 仅占 **0.5%**，几乎可忽略。

---

## 四、具体要修改 FlowDock 的哪些代码

### 4.1 FlowDock 代码结构分析

FlowDock 的关键推理代码位于以下文件中（需克隆后确认具体路径）：

| 文件 | 功能 | 修改需求 |
|:---|:---|:---|
| `inference.py` / `sample.py` | 推理入口，调用 ODE 求解器 | **主修改点**：在 ODE 循环中插入投影 |
| `model.py` / `flow_model.py` | 神经网络定义 | **不修改** |
| `data.py` / `dataset.py` | 数据加载、分子图解析 | **小修改**：提取配体邻接矩阵、键长/键角范围 |
| `utils.py` | 工具函数 | **不修改** |
| `openmm_relax.py` | OpenMM 后处理 | **不修改**（但在实验中可以跳过或对比） |

### 4.2 需要新增的文件

| 新文件 | 功能 |
|:---|:---|
| `constraints.py` | 约束定义（C1-C5 的残差函数） |
| `projection.py` | CCFM 投影算子（第一层 + 第二层 BCP） |
| `scheduler.py` | 概率调度器 $\phi(t)$ |

### 4.3 修改清单

**修改 1：数据预处理阶段——提取约束参数**

在数据加载时，从配体 SMILES 中提取：

```python
# 新增函数：从 RDKit mol 对象中提取约束信息
def extract_constraint_params(mol):
    """
    Args:
        mol: RDKit Mol 对象
    Returns:
        bond_pairs: list of (i, j) 共价键对
        bond_bounds: list of (lower, upper) 键长范围 (Å)
        angle_triples: list of (i, j, k) 键角三元组
        angle_bounds: list of (lower, upper) 键角范围 (rad)
        nonbond_pairs: list of (i, j) 非键合原子对 (拓扑距离 >= 3)
        vdw_radii: array of van der Waals radii per atom
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors

    # 1. 提取共价键对及键长范围
    bond_pairs = []
    bond_bounds = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond_pairs.append((i, j))
        # 使用 RDKit 的 GetBondLength 或 CSD 数据
        bt = bond.GetBondType()
        ai = mol.GetAtomWithIdx(i).GetSymbol()
        aj = mol.GetAtomWithIdx(j).GetSymbol()
        low, up = get_bond_length_bounds(ai, aj, bt)  # 查表函数
        bond_bounds.append((low, up))

    # 2. 提取键角三元组及范围
    angle_triples = []
    angle_bounds = []
    for atom in mol.GetAtoms():
        j = atom.GetIdx()
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        for a in range(len(neighbors)):
            for b in range(a+1, len(neighbors)):
                i, k = neighbors[a], neighbors[b]
                angle_triples.append((i, j, k))
                low, up = get_angle_bounds(mol, i, j, k)  # 查表函数
                angle_bounds.append((low, up))

    # 3. 提取非键合原子对（排除 1-2 和 1-3 连接）
    # 通过图最短路径距离判断
    dm = Chem.GetDistanceMatrix(mol)
    nonbond_pairs = []
    N = mol.GetNumAtoms()
    for i in range(N):
        for j in range(i+1, N):
            if dm[i][j] >= 3:
                nonbond_pairs.append((i, j))

    # 4. 范德华半径
    vdw_radii = np.array([
        Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum())
        for atom in mol.GetAtoms()
    ])

    return bond_pairs, bond_bounds, angle_triples, angle_bounds, nonbond_pairs, vdw_radii
```

**修改 2：推理循环——核心修改**

定位 FlowDock 推理代码中的 ODE 循环（通常形如 `for n in range(N_steps):`），在每步 VD-ODE 更新后插入投影：

```python
# ===== 修改位置：FlowDock inference.py 中的 ODE 循环 =====

from physflowdock.projection import ligand_internal_projection, bilateral_clash_projection
from physflowdock.scheduler import PhiScheduler

# 初始化调度器
phi_scheduler = PhiScheduler(t_switch=0.5, n1=0.5, n2=0.1)

# 初始化约束参数（在推理开始前，从配体信息中提取一次）
constraint_params = extract_constraint_params(ligand_mol)
prot_vdw = get_protein_vdw_radii(protein_atoms)

x = x_0  # 初始状态 [x_P, x_L], shape: (N_prot + N_lig, 3)
for n in range(N_steps):
    t = n / N_steps
    s = (n + 1) / N_steps

    # --- FlowDock 原始 VD-ODE 更新（保持不变）---
    v = model(x, t)
    alpha = torch.clamp(((1 - s) / (1 - t)) ** eta, 0.0, 1.0)
    x_tilde = alpha * x + (1 - alpha) * v

    # --- PhysFlowDock 新增：物理约束投影 ---
    phi_t = phi_scheduler(s)  # 使用 s = t_{n+1}

    if phi_t > 0.01:  # 避免在 t ≈ 0 时浪费计算
        # Step 1: 配体内部约束
        x_tilde = ligand_internal_projection(
            x_tilde, phi_t, s,
            constraint_params=constraint_params,
            max_iter=2, lam=1e-6
        )

        # Step 2: 蛋白-配体双向协同投影
        x_tilde = bilateral_clash_projection(
            x_tilde, phi_t, s,
            prot_vdw=prot_vdw,
            lig_vdw=constraint_params['vdw_radii'],
            delta_pl=0.75,
            max_iter=2, lam=1e-6,
            w_P=0.3
        )

    x = x_tilde

# 最终输出 x 即为物理一致的复合物结构
```

**修改 3：不需要修改的部分（明确声明）**

- **模型权重**：不修改。使用 FlowDock 官方预训练权重。
- **训练代码**：不修改。PhysFlowDock 是 training-free 方法。
- **损失函数**：不修改。
- **网络架构**：不修改。

---

## 五、与 CCFM 的实验对比方案

### 5.1 对比矩阵

| 方法 | 骨干模型 | 蛋白表示 | 碰撞投影策略 | 后处理 |
|:---|:---|:---|:---|:---|
| FlexDock | FlexDock | 扭转角 | 无 | 无 |
| CCFM + FlexDock | FlexDock | 扭转角 | 单向（仅动配体） | 无 |
| FlowDock（原始） | FlowDock | 全原子笛卡尔 | 无 | OpenMM |
| FlowDock（无后处理） | FlowDock | 全原子笛卡尔 | 无 | 无 |
| **PhysFlowDock（$w_P=0$）** | FlowDock | 全原子笛卡尔 | 单向（仅动配体） | 无 |
| **PhysFlowDock（$w_P=0.3$）** | FlowDock | 全原子笛卡尔 | **双向（BCP）** | 无 |
| **PhysFlowDock + 轻量 OpenMM** | FlowDock | 全原子笛卡尔 | **双向（BCP）** | 少量迭代 |

其中 **PhysFlowDock ($w_P=0$)** 是关键消融：它等价于"CCFM 思想移植到 FlowDock 的 VD-ODE"，**去掉了双向投影的创新**。与 PhysFlowDock ($w_P=0.3$) 的对比可以直接量化双向协同投影的增益。

### 5.2 核心指标

| 指标 | 说明 |
|:---|:---|
| RMSD ≤ 2Å (%) | 对接精度 |
| PB-Valid (%) | PoseBusters 物理合法性 |
| PB-Valid ∩ RMSD ≤ 2Å (%) | 联合成功率（最核心指标） |
| 推理时间 (s) | 总时间（含投影，不含 OpenMM） |
| TPPI | 轨迹物理合理性指数 |

### 5.3 关键消融实验

| 消融 | 对比组 | 验证目标 |
|:---|:---|:---|
| **BCP 有效性** | $w_P=0$ vs $w_P=0.3$ vs $w_P=0.5$ vs $w_P=1.0$ | 双向投影的最优权重 |
| **约束层次** | 仅 C1-C3 vs 仅 C4-C5 vs 全部 | 每层约束的边际贡献 |
| **调度器** | $n_1=n_2=0.5$ / $n_1=n_2=0.1$ / 两阶段 | 调度策略的影响 |
| **投影迭代次数** | 1 / 2 / 3 / 5 次 Gauss-Newton | 收敛所需的最小迭代 |
| **OpenMM 协同** | PhysFlowDock + 0/10/50/100 步 OpenMM | 投影后是否还需要后处理 |

---

## 六、论文 Storyline

### 标题候选

**PhysFlowDock: Bilateral Cooperative Projection for Physically Consistent Protein–Ligand Docking via Chance-Constrained Flow Matching**

### 核心卖点

1. **Bilateral Cooperative Projection（双向协同投影）**：首次在生成式对接中实现蛋白-配体双向避让的物理约束投影，模拟真实的诱导契合效应。与 CCFM 的"单向推配体"形成鲜明对比。

2. **VD-ODE 适配**：首次将 CCFM 理论从标准 Euler FM 推广到 Variance Diminishing ODE 求解器，并设计了适配 VD-ODE 的经验概率调度器。

3. **消除 OpenMM 依赖**：在 FlowDock（当前最强的全柔性对接模型）上实现了免后处理的物理一致性生成，开销仅增加 ~0.2 秒。

4. **Physical Navigation（物理导航）**：配体不再"幽灵穿墙"，而是在蛋白表面寻找空隙"滑入"口袋。蛋白残基也会"让位"，展现诱导契合的动态过程。

---

## 七、工程实施路线图

### Week 1：环境搭建 + Baseline 数据收集

- 克隆 FlowDock 仓库，配置环境
- 在 PoseBusters Benchmark (308) 上运行 FlowDock **不加 OpenMM**，收集 Baseline 1
- 在 PoseBusters Benchmark 上运行 FlowDock **加 OpenMM**，收集 Baseline 2
- **定位 VD-ODE 循环的精确代码位置**

### Week 2：实现约束定义 + 第一层投影

- 实现 `constraints.py`：C1-C3 残差函数
- 实现 `projection.py`：第一层投影（配体内部）
- 实现 `scheduler.py`：概率调度器
- 选 10 个复合物快速验证：检查键长/键角是否改善

### Week 3：实现第二层投影（BCP）

- 实现 `projection.py`：第二层双向协同投影
- 调试加权 Gauss-Newton 更新
- 选 10 个复合物验证：检查碰撞是否消除、蛋白是否合理移动

### Week 4：全量实验 + 消融

- PoseBusters Benchmark 全量运行（四组实验）
- 运行关键消融实验（$w_P$ 消融最优先）
- 收集推理时间数据

### Week 5-6：可视化 + 论文撰写

- 制作轨迹对比动画
- 计算 TPPI
- 撰写论文

---

## 八、预期成果

| 指标 | FlowDock (无后处理) | FlowDock + OpenMM | CCFM + FlexDock (10 steps) | **PhysFlowDock (预期)** |
|:---|:---|:---|:---|:---|
| PB-Valid (%) | ~41% | ~51% | ~78.8% | **≥ 80%** |
| PB-Valid ∩ L-RMSD < 2Å (%) | — | ~51% | ~30.8% | **≥ 45%** |
| 推理总时间 (s) | ~39 | ~39 + 数十秒 | ~0.79 | **~40** |
| 需要 OpenMM？ | 必须 | 是 | 否 | **否** |
| 蛋白碰撞投影 | 无 | 后处理 | 单向（仅动配体） | **双向（BCP）** |

注意：CCFM+FlexDock 的推理极快（<1s），是因为 FlexDock 使用低维参数空间且网络较小。FlowDock 的推理时间~39s 来自其大型 SE(3)-等变 Transformer + ESMFold。PhysFlowDock 的额外开销（~0.2s）相对于 39s 可忽略，但总时间不可能降到 1s 以内（这是骨干模型的固有成本）。
