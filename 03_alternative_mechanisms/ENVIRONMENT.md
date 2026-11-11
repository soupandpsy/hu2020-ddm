# M0 / M1 运行环境

分析使用 Linux 下的 dockerHDDM：HDDM 1.0.1RC、Kabuki 0.6.5RC4、ArviZ 0.15.1 和原生 PyMC 2 采样器。Notebook 文件可以位于 Windows；执行环境由所选内核决定。

原运行镜像为 `hu2020-part3:pyddm090-tail`，基于 `hcp4715/hddm:1.0.1`。新增的 PyDDM 动态求解组件不参与 M0/M1 的原生密度、拟合和 PPC。这两份 Notebook 可使用原始 HDDM 镜像运行；不需要 M3 的 JAX/HSSM 环境。

基础镜像版本固定为：

```text
hcp4715/hddm@sha256:5b30f195d3d41905f5461dbe9c786ea2ab531b4355d7464cb0fbfb5d6adfd238
```

在仓库根目录的 PowerShell 中启动 Jupyter：

```powershell
$projectRoot = (Get-Location).Path
docker run --rm -it -p 127.0.0.1:8888:8888 -v "${projectRoot}:/repo" -w /repo hcp4715/hddm@sha256:5b30f195d3d41905f5461dbe9c786ea2ab531b4355d7464cb0fbfb5d6adfd238 jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```

打开终端给出的带 token 链接，进入 `03_alternative_mechanisms/Scripts/notebooks/`。Notebook 第一章显示实际 Python 和依赖版本。保留默认身份验证。

## 执行与结果读取

所有 `RUN_*` 计算开关默认为 `False`。不需要重新运行即可阅读保存在 Notebook 中的结果。新克隆没有原始链，直接 Run All 不会自动采样，但会将部分结果输出更新为“缺少完整链后验”；查看发布结果时应保留原保存输出。

重新分析时，两个 Notebook 分别执行 M0 和 M1，不依赖彼此的内存状态。预算为每模型四链，每链 10,000 次、burn-in 5,000、thinning 1。`N_WORKERS=1` 控制单页并发；全机并发上限为 16。M0 种子为 2026110719，M1 为 2026110791。

`RUN_FIT` 控制拟合；`RUN_POSTPROCESS` 保存诊断和参数概要；`RUN_DIC`、`RUN_PPC` 和 `RUN_PPC_DENSITY` 分别计算 DIC、行为统计量 PPC 和 RT 密度图。`RUN_PRIOR` 仅抽取联合参数先验，本次报告没有执行该 2,000 次先验抽样。PPC 总量为 500×完整链数，四链合计 2,000 个复制数据集；RT 密度图重放其中前 100 组后验抽样位置及随机种子。

完整链、DIC 缓存和 PPC 记录位于：

```text
Output/tables/local_bulk/chain_comparison_20261109/M0/
Output/01_starting_point_z/local_bulk/chain_comparison_20261109/M1/
```

已有完整链按模型、输入、预算和源代码指纹核验后复用。原始 `.db`、`.hddm`、`.nc` 及 `local_bulk/` 未发布；恢复本地档案时需保持对应目录结构。仅有汇总 CSV 不能恢复链轨迹或重新计算 PPC。重新拟合可使用新的 `RUN_ID` 保留原结果。

## 数据与代码版本

两个模型共同使用 19,580 条固定训练试次，哈希为 `2dfcb00ea2bc985b516d40614a8af2cde3cbc004f96d67fd2f7ef2ff13456951`。四个共享 Python 文件的 SHA256 记录在各模型 provenance 中，LF 换行由 `.gitattributes` 固定，防止跨系统换行变化被误判为模型修改。

共享文件包含其他机制的辅助接口，但此发布版本的 M0/M1 入口不调用动态求解器。模型、先验和拟合结果保持原值；既有 RT 筛选与 Deadline 限制见两份 Notebook。
