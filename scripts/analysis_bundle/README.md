# design031726 — Analysis scripts bundle

本目录为 **HBV design031726** 项目下与「合并打分 / KIH packing / 相关性图 / Top 筛选」相关的脚本副本，便于拷贝到其它机器或归档。

## 依赖

```bash
module load python   # 按集群习惯
pip install --user -r requirements.txt
```

## 脚本说明

| 脚本 | 作用 |
|------|------|
| `analyze_merged_scores.py` | 合并 MPNN CSV + Rosetta scorefile，算 composite，出多种图 |
| `analyze_kih_packing.py` | 对 Rosetta 输出 PDB 做简化 KIH（knob/hole）分析，批量 CSV |
| `analyze_rmsd.py` | AF3 vs Rosetta relaxed 的 backbone RMSD |
| `plot_knob_correlation_heatmap.py` | 合并 merged + KIH CSV，画相关性热图（可小图 subset） |
| `plot_knob_histograms.py` | KIH knob 相关指标直方图 |
| `copy_knob_gt_threshold.py` | 按 `knob_score` 阈值复制 PDB 到新目录 |
| `select_top20_kih_combo.py` | 按 sc / knob / i_ptm / rmsd / holes 等综合分选 TopN 并复制 PDB |

## 路径约定（示例）

- Rosetta relaxed PDB 目录：`.../outputs_4/rosetta_fastrelax/`
- MPNN+Rosetta 合并表：`.../outputs_4/merged_scores.csv`
- KIH 结果：`.../outputs_4/kih_packing_results_opt.csv`

运行前请把上述路径改成你机器上的绝对路径。

## 打包整个目录为压缩包

在 `design031726` 上一级执行：

```bash
tar czvf design031726_analysis_bundle.tgz design031726/analysis_bundle
```

或在 `design031726` 内：

```bash
tar czvf analysis_bundle.tgz analysis_bundle
```
