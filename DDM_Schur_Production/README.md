# Dynamic Schur Production

这是从研究仓库中整理出的独立生产工程，只保留目前在 Package15
chiplet 原型上最快且已验证精度的组合：

- 区域内部：`M1 Local Block-Arnoldi`；
- 构造激励：`global-fom` operator-derived traces；
- 界面：保留完整物理界面自由度，不做 POD 或谱端口降维；
- 在线求解：`augmented-direct`，固定时间步下复用一次 MKL PARDISO
  因式分解；
- 时间推进：Backward Euler + native reduced history；
- 校验策略：不做残差回退，不输出完整场；可选 summary-only 单体 SIPG
  校验，且 FOM 不参与 ROM 构造。

研究仓库没有被修改。本目录有自己的 CMake、CLI、配置、网格、测试和文档。

## 快速开始

环境要求：Visual Studio 2019 或更新版本、CMake 3.20+、Intel oneAPI MKL。

```powershell
cd .\DDM_Schur_Production
.\scripts\build_release.ps1 -Jobs 16
.\scripts\run_package15.ps1 -Profile smoke -Stage Both -Steps 2 -Force
ctest --test-dir build -C Release --output-on-failure
```

`Cold` 会生成 descriptor/reference/local-ROM 缓存；`Warm` 要求缓存已经存在，
并拒绝静默重建。输出目录必须为空，避免把不同配置的结果混在一起。

大模型单步示例：

```powershell
.\scripts\run_package15.ps1 `
  -Profile large -Stage Cold -Steps 1 -Threads 16 -LocalMklThreads 1 -Force
```

仅在最终精度审计时增加单体 SIPG summary-only 对照：

```powershell
.\scripts\run_package15.ps1 `
  -Profile large -Stage Warm -Steps 1 -ValidateWithFom -Force
```

## 自带模型

| Profile | 子域 | 界面 | 独立 P2 DOFs | 用途 |
|---|---:|---:|---:|---|
| `smoke` | 15 | 18 | 21,507 | 数秒级回归测试 |
| `medium` | 15 | 18 | 92,355 | 日常开发验证 |
| `large` | 15 | 18 | 503,209 | 工作站性能评估 |

15 个子域由 8 个同构 HBM、1 个 Logic、4 个 Interposer tile、1 个
Substrate 和 1 个 Baseplate 组成。HBM 使用共享 ROM template，但每个实例
保留独立全局映射；这对一个此前没有计算过、但包含重复几何/材料组件的新模型
同样是合法的结构复用。

## 已验证基线

整理前同一数值路径在 Package15 large 上得到：503,209 全局 DOFs、214,700
界面 DOFs、总局部秩 540；冷启动由 27.3361 s 降至 22.1535 s，单步约
0.3127 s。summary-only FOM 对照为相对 L2 误差 `2.60e-15`、最大温度误差
`1.48e-12 K`、完整残差 `6.91e-8`。这些是已保存的历史基线，不冒充本次
源码整理后的重新跑分；本工程当前已完成独立 Release 和 smoke 冷/热验证。

## 目录

- `src/main.cpp`：受限 CLI 和固定生产参数；
- `src/mor/transient/local_dynamic_schur.cpp`：descriptor、缓存、Arnoldi、
  时间推进和汇总；
- `src/mor/local/local_reduced_schur.cpp`：仅保留 augmented-direct；
- `src/internal/removed_research_methods.cpp`：被删除研究方法的不可达保护边界；
- `configs/`：三个 Package15 配置；
- `data/generated/package15/`：已生成的可运行网格；
- `tools/generate_package15_mesh.py`：确定性网格生成器；
- `tests/run_production_smoke.ps1`：冷缓存、热缓存与算法契约回归；
- `docs/ALGORITHM_AND_CODE.md`：算法、模块、数据流和流程图；
- `docs/REMOVED_METHODS.md`：删除范围与理由。

主要运行结果是 `local_dynamic_schur_summary.csv`、
`local_dynamic_schur_timing.csv`、`augmented_direct_summary.csv` 和
`production_run.csv`。生产入口不会写 `local_dynamic_schur_final_temperature.csv`。
