# SynComDesign 中文使用手册

本手册面向没有 MATLAB、COBRA Toolbox、代谢模型或 FBA 基础的用户。按照下面步骤操作，可以完成模型准备、培养基设置、群落组合分析和结果查看。

## 1. 软件功能简介

SynComDesign 是一个用于合成微生物群落设计和功能评价的 MATLAB 工具。它可以读取多个菌株的基因组尺度代谢模型，自动枚举菌株组合，构建共享环境群落模型，并预测生长和反硝化相关通量。

当前版本支持：

- CarveMe 生成的 COBRA/SBML 模型。
- BiGG 风格的反应和代谢物 ID。
- 单菌株和多菌株组合枚举。
- 群落总生物量预测。
- 指定目标菌生物量预测。
- 等比例或固定比例群落组成约束。
- 生长优先再优化 N2O 消耗的二阶段目标。
- 硝酸盐、亚硝酸盐、一氧化氮、氧化亚氮和氮气交换通量输出。
- TSV 和 CSV 格式结果表。

## 2. 文件夹结构

GitHub 上传版推荐结构如下：

```text
SynComDesign/
  README.md
  .gitignore
  config/
    syncomdesign_config.yml
    biomass_reactions.tsv
    metabolite_aliases.tsv
  media/
    medium.tsv
  models/
    005.xml
    016.xml
    020.xml
  scripts/
    runSynComDesign.m
    GetAllCombination.m
    src/  
  USER_MANUAL_CN.md
```

各文件夹含义：

- `models/`：放菌株代谢模型。
- `media/`：放培养基文件。
- `config/`：放主配置文件和映射表。
- `scripts/`：放 MATLAB 程序。
- `USER_MANUAL_CN.md`：中文使用说明。

## 3. 需要安装的软件

运行 SynComDesign 需要：

1. MATLAB。
2. COBRA Toolbox。参考"https://opencobra.github.io/cobratoolbox/stable/installation.html".
3. 一个 LP 求解器，例如 Gurobi 或 GLPK。

推荐使用 Gurobi。

## 4. 第一次运行

打开 MATLAB，进入项目目录：

```matlab
cd('path\to\SynComDesign')
addpath(genpath(pwd))
```

把 `path\to\SynComDesign` 换成你自己电脑上的实际路径。

初始化 COBRA Toolbox 和求解器：

```matlab
initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
```

运行 SynComDesign：

```matlab
results = runSynComDesign('config/syncomdesign_config.yml');
```

如果成功，会看到类似输出：

```text
Models detected: 3
Models validated: 3
Biomass reactions detected: 3
Total combinations: 7
Feasible combinations: 7
Failed combinations: 0
Output directory: ...\results
```

## 5. 主配置文件

主配置文件是：

```text
config/syncomdesign_config.yml
```

常用配置如下：

```yaml
project:
  name: SynComDesign
  output_dir: results

models:
  directory: models
  file_pattern: "*.xml"
  biomass_reactions_file: config/biomass_reactions.tsv
  metabolite_aliases_file: config/metabolite_aliases.tsv

medium:
  file: media/medium.tsv
  condition: anaerobic
  close_unlisted_uptakes: true

objective:
  scenario_id: 1
  growth_fraction: 0.9
  target_strain: null
  biomass_weights: equal

solver:
  name: gurobi
```

## 6. 如何设置输出文件夹

结果保存位置由这里控制：

```yaml
project:
  output_dir: results
```

如果不想覆盖旧结果，可以改成：

```yaml
output_dir: results_ID1_growth
```

每次换一个 `output_dir`，结果就会保存到不同文件夹。

## 7. 如何放置模型

模型放在：

```text
models/
```

例如：

```text
models/005.xml
models/016.xml
models/020.xml
```

配置文件中对应：

```yaml
models:
  directory: models
  file_pattern: "*.xml"
```

如果你换成自己的模型，建议模型文件名直接使用菌株名，例如：

```text
models/strainA.xml
models/strainB.xml
models/strainC.xml
```

模型文件名去掉 `.xml` 后，就是结果表里的菌株名。

## 8. biomass_reactions.tsv

文件位置：

```text
config/biomass_reactions.tsv
```

它告诉程序每个模型的 biomass 反应是哪一个。

格式示例：

```text
strain	biomass_rxn
005	Growth
016	Growth
020	Growth
```

第一列 `strain` 必须和模型文件名一致。比如模型是：

```text
models/005.xml
```

那么这里写：

```text
005
```

第二列 `biomass_rxn` 必须是模型中真实存在的 biomass reaction ID。

如果 biomass reaction 写错，结果可能全失败，或者生物量为 0。

## 9. 如何设置培养基

培养基文件放在：

```text
media/medium.tsv
```

配置文件中对应：

```yaml
medium:
  file: media/medium.tsv
```

COBRA 通常使用负数表示摄取。例如：

```text
EX_no3_e	-10	1000
```

表示允许硝酸盐摄取，最大摄取速率为 10。

如果你想测试新培养基，建议不要覆盖原文件，而是复制一份：

```text
media/medium_low_nitrate.tsv
```

然后修改配置：

```yaml
medium:
  file: media/medium_low_nitrate.tsv
```

## 10. metabolite_aliases.tsv

文件位置：

```text
config/metabolite_aliases.tsv
```

它告诉程序哪些模型 ID 对应硝酸盐、亚硝酸盐、NO、N2O 和 N2。

示例：

```text
canonical_id	alias	category
n2o	n2o_e	nitrous_oxide
n2o	EX_n2o_e	nitrous_oxide
```

如果你的模型使用不同的 N2O 反应或代谢物 ID，需要把它加到这个表里。

这个文件不建议删除。删除后程序可能无法正确识别反硝化通量。

## 11. ID1 到 ID5 怎么选

目标函数由 `scenario_id` 控制：

```yaml
objective:
  scenario_id: 1
```

### ID1：最大化群落总生物量

```yaml
objective:
  scenario_id: 1
```

用途：看哪个组合整体长得最好。

### ID2：最大化指定目标菌生物量

```yaml
objective:
  scenario_id: 2
  target_strain: E10
```

用途：看哪些组合最有利于目标菌生长。

注意：`target_strain` 必须和模型文件名一致。如果模型是 `E10.xml`，这里写 `E10`。

如果目标菌生物量仍然是 0，说明该菌在当前培养基和模型约束下可能不能生长。

### ID3：等比例群落组成

```yaml
objective:
  scenario_id: 3
```

用途：要求组合中菌株 biomass 通量保持等比例。

### ID4：固定比例群落组成

```yaml
objective:
  scenario_id: 4
```

如果没有设置具体比例，当前程序会按默认等比例处理。

### ID5：兼顾生长和 N2O 消耗

```yaml
objective:
  scenario_id: 5
  growth_fraction: 0.9
```

含义：

1. 先最大化生长。
2. 至少保留 90% 最大生长量。
3. 在这个前提下尽量提高 N2O 摄取。

如果当前模型没有 N2O 摄取能力，ID5 的结果可能和 ID1 一样。

## 12. 结果文件

运行结束后会生成结果文件夹，例如：

```text
results/
```

主要文件：

```text
community_summary.tsv
community_summary.csv
single_strain_results.tsv
model_validation.tsv
reaction_mapping.tsv
metabolite_mapping.tsv
failed_combinations.tsv
run.log
```

最重要的是：

```text
community_summary.tsv
```

常用列：

- `combination_id`：组合名称。
- `community_size`：组合中菌株数量。
- `strain_names`：组合包含的菌株。
- `feasible`：是否可行，1 表示可行。
- `objective_mode`：使用的目标模式。
- `total_biomass`：总生物量。
- `strain_biomass_*`：每个菌株的生物量。
- `nitrate_uptake`：硝酸盐摄取。
- `nitrite_uptake`：亚硝酸盐摄取。
- `nitrite_secretion`：亚硝酸盐分泌。
- `no_uptake`：一氧化氮摄取。
- `no_secretion`：一氧化氮分泌。
- `n2o_uptake`：氧化亚氮摄取。
- `n2o_secretion`：氧化亚氮分泌。
- `n2o_net_flux`：N2O 净通量。
- `n2_secretion`：氮气产生。

## 13. 常见问题

### 13.1 报错：Default LP solver not selected

说明没有设置 COBRA 求解器。

在 MATLAB 里运行：

```matlab
initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
```

然后重新运行。

### 13.2 所有组合都失败

打开：

```text
results/failed_combinations.tsv
```

查看 `error_message`。

常见原因：

- 求解器没有设置。
- biomass reaction 写错。
- 模型文件路径不对。
- 培养基文件路径不对。
- 模型中缺少关键 exchange reaction。

### 13.3 目标菌生物量为 0

可能原因：

- 目标菌在当前培养基下不能生长。
- 目标菌 biomass reaction 写错。
- 培养基缺少必要营养物。
- 模型本身缺少生长所需路径。

建议先用 ID1 检查该菌单独是否能生长。

### 13.4 N2O uptake 全是 0

可能原因：

- 模型没有 N2O exchange reaction。
- `metabolite_aliases.tsv` 没有写对 N2O ID。
- 培养基没有允许 N2O 摄取。
- 模型没有利用 N2O 的代谢路径。

## 14. 推荐使用流程

每次正式分析建议按这个顺序：

1. 检查 `models/` 中模型是否正确。
2. 检查 `config/biomass_reactions.tsv`。
3. 检查 `media/medium.tsv`。
4. 检查 `config/metabolite_aliases.tsv`。
5. 设置 `scenario_id`。
6. 设置 `output_dir`。
7. 运行 `runSynComDesign`。
8. 查看 `community_summary.tsv`。
9. 查看 `failed_combinations.tsv`。

## 15. 引用

当前项目参考并修改了原 SuperCC 工作流。公开发布前，请保留原论文和 GitHub 仓库引用，并确认相关代码、模型和第三方依赖的许可要求。

原始项目：

```text
https://github.com/ruanzhepu/superCC
```

原始论文：

```text
Ruan, Z., Chen, K., Cao, W. et al. Engineering natural microbiomes toward
enhanced bioremediation by microbiome modeling. Nature Communications 15, 4694
(2024). https://doi.org/10.1038/s41467-024-49098-z
```
