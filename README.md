# SynComDesign

**Constraint-Based Design and Functional Evaluation of Synthetic Microbial Communities**

SynComDesign is a MATLAB-based workflow for designing and evaluating synthetic microbial communities with genome-scale metabolic models. It reads strain models, applies a user-defined medium, enumerates strain combinations, builds shared-environment community models, and reports growth and denitrification-related exchange fluxes.

This README is written as a beginner-friendly guide. You do not need prior experience with SynComDesign, but you do need MATLAB, COBRA Toolbox, and a linear programming solver.

## 1. What SynComDesign Does

SynComDesign can be used to answer questions such as:

- Which strain combination has the highest predicted total biomass?
- Which community best supports a target strain?
- Which combinations remain feasible under equal or fixed composition constraints?
- Which combinations consume nitrate or nitrous oxide?
- Which combinations produce dinitrogen?

Current capabilities include:

- CarveMe-generated COBRA/SBML model support.
- BiGG-style reaction and metabolite identifier support.
- Enumeration of single-strain and multi-strain communities.
- Shared extracellular environment construction.
- Total community biomass prediction.
- Target-strain biomass optimization.
- Equal-composition and fixed-composition community constraints.
- Growth-first, N<sub>2</sub>O-consumption-second optimization.
- Nitrate, nitrite, nitric oxide, nitrous oxide, and dinitrogen exchange reporting.
- TSV and CSV output tables.

## 2. Repository Structure

The GitHub upload package is organized as:

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

Folder meanings:

- `models/`: strain metabolic models.
- `media/`: medium files.
- `config/`: main configuration and mapping files.
- `scripts/`: MATLAB entry point and source code.
- `docs/`: additional documentation.

## 3. Software Requirements

You need:

1. MATLAB.
2. COBRA Toolbox."https://opencobra.github.io/cobratoolbox/stable/installation.html"
3. A COBRA-compatible LP solver, such as Gurobi or GLPK.

Gurobi is recommended when available. GLPK can be used, but it may be less stable for some genome-scale models.

CarveMe and BiGG are not bundled with this repository. SynComDesign is compatible with CarveMe-generated models and BiGG-style identifiers.

## 4. First Run

Open MATLAB and move into the downloaded SynComDesign folder:

```matlab
cd('path\to\SynComDesign')
addpath(genpath(pwd))
```

Replace `path\to\SynComDesign` with the actual path on your computer.

Initialize COBRA Toolbox and select a solver:

```matlab
initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
```

Run SynComDesign:

```matlab
results = runSynComDesign('config/syncomdesign_config.yml');
```

If the run succeeds, MATLAB will print a summary similar to:

```text
Models detected: 3
Models validated: 3
Biomass reactions detected: 3
Total combinations: 7
Feasible combinations: 7
Failed combinations: 0
Output directory: ...\results
```

## 5. Main Configuration File

The main configuration file is:

```text
config/syncomdesign_config.yml
```

Important sections:

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

## 6. Output Folder

The output folder is controlled by:

```yaml
project:
  output_dir: results
```

To avoid overwriting previous runs, use a different output folder for each analysis:

```yaml
output_dir: results_ID1_growth
```

or:

```yaml
output_dir: results_ID2_E10
```

Each run writes result tables into the selected folder.

## 7. Input Models

Put strain models in:

```text
models/
```

Example:

```text
models/005.xml
models/016.xml
models/020.xml
```

The configuration reads XML models by default:

```yaml
models:
  directory: models
  file_pattern: "*.xml"
```

If you use your own models, name the files clearly:

```text
models/strainA.xml
models/strainB.xml
models/strainC.xml
```

The file name without `.xml` becomes the strain name used in output tables.

## 8. Biomass Reaction File

The biomass reaction mapping file is:

```text
config/biomass_reactions.tsv
```

It tells SynComDesign which biomass reaction belongs to each strain.

Example:

```text
strain	biomass_rxn
005	Growth
016	Growth
020	Growth
```

The first column, `strain`, must match the model file name without the extension. For `models/005.xml`, the strain name is `005`.

The second column, `biomass_rxn`, must be a real reaction ID in the model.

Do not delete this file unless you are certain automatic biomass detection works for all models. If the biomass reaction is wrong, runs may fail or report zero growth.

## 9. Medium File

The medium file is:

```text
media/medium.tsv
```

The configuration points to it here:

```yaml
medium:
  file: media/medium.tsv
```

COBRA convention usually uses negative lower bounds for uptake. For example:

```text
EX_no3_e	-10	1000
```

This means nitrate uptake is allowed with a maximum uptake rate of 10.

To test a new medium, copy the medium file instead of overwriting it:

```text
media/medium_low_nitrate.tsv
```

Then update the configuration:

```yaml
medium:
  file: media/medium_low_nitrate.tsv
```

## 10. Metabolite Alias File

The metabolite alias file is:

```text
config/metabolite_aliases.tsv
```

It tells SynComDesign which model IDs correspond to nitrate, nitrite, NO, N<sub>2</sub>O, and N<sub>2<sub>.

Example:

```text
canonical_id	alias	category
n2o	n2o_e	nitrous_oxide
n2o	EX_n2o_e	nitrous_oxide
```

If your model uses a different N<sub>2</sub>O exchange reaction or metabolite ID, add it to this table.

This file should normally be kept. If it is removed, denitrification flux extraction may become incomplete or incorrect.

## 11. Objective Modes: ID1 to ID5

The objective mode is controlled by:

```yaml
objective:
  scenario_id: 1
```

### ID1: Maximize Total Community Biomass

```yaml
objective:
  scenario_id: 1
```

Use this mode to find the combination with the highest total predicted growth.

### ID2: Maximize Target-Strain Biomass

```yaml
objective:
  scenario_id: 2
  target_strain: 005
```

Use this mode when you care about one specific strain. The target strain name must match the model file name. For `E10.xml`, use:

```yaml
target_strain: 005
```

SynComDesign automatically evaluates only combinations that contain the target strain.

If the target-strain biomass is still zero, the target strain may not be able to grow under the current model and medium constraints.

### ID3: Equal Community Composition

```yaml
objective:
  scenario_id: 3
```

This mode constrains strain biomass fluxes to an equal ratio within each community.

### ID4: Fixed Community Composition

```yaml
objective:
  scenario_id: 4
```

This mode is used for fixed biomass-ratio constraints. If no custom ratio is provided, the current implementation falls back to equal ratios.

### ID5: Growth-Then-N<sub>2</sub>O-Consumption

```yaml
objective:
  scenario_id: 5
  growth_fraction: 0.9
```

This mode runs a two-step optimization:

1. Maximize total growth.
2. Keep at least 90% of the maximum growth.
3. Within that growth constraint, maximize N<sub>2</sub>O uptake.

If the models cannot consume N<sub>2</sub>O under the current medium, ID5 may give the same result as ID1.

## 12. Output Files

After a run, SynComDesign creates an output folder such as:

```text
results/
```

Typical files:

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

The most important file is:

```text
community_summary.tsv
```

Important columns:

- `combination_id`: strain combination name.
- `community_size`: number of strains in the community.
- `strain_names`: strains included in the community.
- `feasible`: whether the optimization was feasible; `1` means feasible.
- `objective_mode`: objective mode used for the run.
- `total_biomass`: total community biomass.
- `strain_biomass_*`: biomass flux for each strain.
- `nitrate_uptake`: nitrate uptake.
- `nitrite_uptake`: nitrite uptake.
- `nitrite_secretion`: nitrite secretion.
- `no_uptake`: nitric oxide uptake.
- `no_secretion`: nitric oxide secretion.
- `n2o_uptake`: nitrous oxide uptake.
- `n2o_secretion`: nitrous oxide secretion.
- `n2o_net_flux`: raw N<sub>2</sub>O exchange direction.
- `n2_secretion`: dinitrogen secretion.

## 13. Flux Sign Convention

Output columns use positive biological rates:

- Uptake is reported as `max(0, -exchangeFlux)`.
- Secretion is reported as `max(0, exchangeFlux)`.
- Net flux keeps the raw COBRA exchange direction.

For example, if the COBRA exchange flux through an N<sub>2</sub>O exchange reaction is `-3`, SynComDesign reports:

```text
n2o_uptake = 3
```

## 14. Common Problems

### 14.1 Error: Default LP Solver Not Selected

Run:

```matlab
initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
```

Then run SynComDesign again.

### 14.2 All Combinations Failed

Open:

```text
results/failed_combinations.tsv
```

Check the `error_message` column.

Common causes:

- The LP solver was not selected.
- A biomass reaction ID is wrong.
- The model path is wrong.
- The medium path is wrong.
- A model is missing key exchange reactions.

### 14.3 Target-Strain Biomass Is Zero

Possible causes:

- The target strain cannot grow in the current medium.
- The target strain biomass reaction is incorrect.
- The medium lacks required nutrients.
- The model lacks a pathway required for growth.

Try ID1 first and check whether the target strain can grow alone.

### 14.4 N<sub>2</sub>O Uptake Is Always Zero

Possible causes:

- The model has no N<sub>2</sub>O exchange reaction.
- `metabolite_aliases.tsv` does not include the correct N<sub>2</sub>O ID.
- The medium does not allow N<sub>2</sub>O uptake.
- The model has no pathway for using N<sub>2</sub>O.

## 15. Recommended Workflow

For each analysis:

1. Check the model files in `models/`.
2. Check `config/biomass_reactions.tsv`.
3. Check `media/medium.tsv`.
4. Check `config/metabolite_aliases.tsv`.
5. Set `scenario_id`.
6. Set `output_dir`.
7. Run `runSynComDesign`.
8. Inspect `community_summary.tsv`.
9. Inspect `failed_combinations.tsv`.

## 16. GitHub Upload Notes

Before uploading to GitHub:

- Do not upload private or restricted models.
- Do not upload generated `results/` folders unless you intend to share example results.
- Keep `config/biomass_reactions.tsv` and `config/metabolite_aliases.tsv`.
- Keep `media/medium.tsv`.
- Keep `scripts/src/`; it contains the core implementation used by `runSynComDesign.m`.

## 17. Notice

SynComDesign was created by modifying and extending an existing implementation of the SuperCC workflow.

The software has been adapted to support CarveMe-generated metabolic models, BiGG identifiers, exhaustive microbial community combination enumeration, shared-environment community models, growth and functional optimization, flux variability analysis, and denitrification-related exchange predictions.

The original SuperCC publication and GitHub repository should be cited when using this software.

SynComDesign is an independently maintained adaptation and is not an official release of the original SuperCC project.

## 18. Citation

If you use SynComDesign, please cite the SynComDesign software repository and the publication describing the original SuperCC framework that informed the development of this project.

Original SuperCC publication:

```text
Ruan, Z., Chen, K., Cao, W. et al. Engineering natural microbiomes toward
enhanced bioremediation by microbiome modeling. Nature Communications 15, 4694
(2024). https://doi.org/10.1038/s41467-024-49098-z
```

Original SuperCC repository:

```text
https://github.com/ruanzhepu/superCC
```

The current citation information for SynComDesign will be added after the software repository or software paper is formally released.
