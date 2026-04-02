# 6QZH 打分结果文件（与 `pdbs/` 分离）

总览与图表说明见 **[README.md](README.md)**（含 **KIH_score** = TSV 列 `packing_finalscore`、**三螺旋 ver3** 含义）。

所有 **packing / Rosetta 合并** 的表格与日志默认放在 **本仓库根目录**（`6QZH/`，与 `scripts/` 同级）；结构文件仍在 **`pdbs/`**（大目录，默认不纳入 Git）。

| 文件 | 说明 |
|------|------|
| `packing_ver3_batch_summary.tsv` | 批量 `pack_analysis_muOR_ver3`：`pdb_path`, `finalscore`, `n_events`, `error` |
| `packing_ver3_batch_errors.log` | 批量运行时的 Python 异常（若有） |
| `packing_scores_ver3.tsv` | 由上表整理的排序版（按 finalscore 降序） |
| `rosetta_packing_merged_scores.tsv` | Rosetta 末尾指标 + packing 列合并 |
| `batch_pack_ver3_stdout.log` | 某次批量运行的终端日志（若存在） |

## 脚本（均在 `scripts/`）

| 脚本 | 作用 |
|------|------|
| `scripts/batch_pack_ver3_all_pdbs.py` | 扫描 `pdbs/**/*.pdb`，结果写入 **仓库根目录**（可用 `MUSOR_SCORE_OUT` 改） |
| `scripts/generate_packing_scores_from_batch_summary.py` | 生成 `packing_scores_ver3.tsv` |
| `scripts/merge_rosetta_and_packing_scores.py` | 读 `pdbs/` 下 PDB 末尾 Rosetta 行 + `packing_scores_ver3.tsv` → `rosetta_packing_merged_scores.tsv` |

示例（在 **6QZH 仓库根目录** 执行）：

```bash
cd /path/to/6QZH
MUSOR_BATCH_JOBS=8 python3 scripts/batch_pack_ver3_all_pdbs.py
python3 scripts/generate_packing_scores_from_batch_summary.py
python3 scripts/merge_rosetta_and_packing_scores.py
```

## 画图（Rosetta + packing 合并表）

依赖：`pip install -r requirements_plots.txt`（或复用 `scripts/analysis_bundle/requirements.txt`）。

```bash
python3 scripts/plot_rosetta_packing_merged.py
# 输出目录默认: figures_rosetta_packing/
#   correlation_heatmap.png
#   scatter_packing_vs_rosetta.png
#   hist_packing_finalscore.png
#   hist_rosetta_metrics.png   (dG_separated, sc_value, Holes, pack_stat, dSASA_int, delta_unsatHbonds)
```

仅使用 Rosetta 行齐全的构象：`python3 scripts/plot_rosetta_packing_merged.py --only-rosetta-ok`

## Sequence logo（PDB 系综）

依赖：`pip install logomaker biopython`（或见 `requirements_plots.txt`）。

从 `pdbs/**/*.pdb` 读取链 **B**（默认）在 **PDB 编号** 上的氨基酸分布，画 Schneider–Stephens 风格 logo。默认输出到 `figures_rosetta_packing/sequence_logo.png`。若要链 A：`--chain A`。

```bash
python3 scripts/plot_sequence_logo.py
# 同一序列只保留一个结构（推荐，581 条 vs 2386 文件）：
python3 scripts/plot_sequence_logo.py --dedupe-by-sequence
# 只画某段残基（例如 100–250）：
python3 scripts/plot_sequence_logo.py --dedupe-by-sequence --region 100 250 --out figures_rosetta_packing/sequence_logo_100-250.png
```
