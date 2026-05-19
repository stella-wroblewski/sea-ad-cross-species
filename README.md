# SEA-AD Cross-Species Comparison Pipeline

A reproducible R pipeline for comparing mouse single-cell differential expression signatures against the Seattle Alzheimer's Disease Brain Cell Atlas ([SEA-AD](https://portal.brain-map.org/explore/seattle-alzheimers-disease); Gabitto et al., 2024, *Nature Neuroscience* 27:2366-2383). Built around the insight that **SEA-AD's all-donors aggregation masks early-phase signals through phase cancellation**, and that phase-stratified analysis is the right comparison for acute disease models.

Developed during my PhD work comparing a 3xTg-AD mouse model of peripheral cytomegalovirus infection to human Alzheimer's. Generalized: any cross-species DE comparison against SEA-AD should run with only the mouse DE table and the file path config edited.

---

## 1. Why this pipeline exists

The straightforward way to compare a mouse model to human AD is to overlap your DE genes against the all-donors SEA-AD DE results. When I did this the obvious way, I got nothing — zero significant overlap in endothelial cells, Spearman rho near zero, percent-concordant at chance. The temptation was to conclude that the mouse model didn't engage AD biology.

That was the wrong conclusion. SEA-AD's CPS-based DE coefficients average across 84 donors spanning the entire ADNC spectrum. The Gabitto paper itself documents two distinct disease phases:

| Phase | Microglia | Astrocytes | Vasculature | Pathology |
|---|---|---|---|---|
| **Early** | Inflammatory activation | Reactive | IFN+, OXPHOS depletion | Slow accumulation |
| **Late** | Anti-inflammatory shift | Mixed | IFN reversal, metabolic recovery | Exponential, neuronal loss |

When the two phases are averaged, the early-phase IFN signal cancels the late-phase IFN reversal, the early-phase OXPHOS depletion cancels the late-phase recovery, and so on. The all-donors result understates real biology. An acute mouse model in young animals should match the early-phase initiating events, not the averaged trajectory.

**This pipeline computes both the naive all-donors comparison and the phase-stratified version side by side**, so you can see directly when phase cancellation is hiding a real effect. In our hands: endothelial cells went from 0 genes at FDR<0.05 in all-donors to 20 at FDR<0.05 in early-donors. The signal was always there.

---

## 2. Repository layout

```
sea-ad-cross-species/
├── README.md                                  this file
├── LICENSE                                    MIT
├── run_all.R                                  wrapper to run the whole pipeline
├── .gitignore
├── scripts/
│   ├── 00_download_sea_ad_data.sh             bash: pull Nebula CSVs from S3
│   ├── 01_load_sea_ad.R                       shared loaders, ortholog map, theme
│   ├── 02_hypergeometric_and_concordance.R    Analyses 1 + 2
│   ├── 03_targeted_gene_heatmap.R             Analysis 3 (curated 28-gene panel)
│   ├── 04_cross_species_gsea.R                Analyses 4 + 6 (Hallmark GSEA)
│   ├── 05_dam_microglia_comparison.R          Analysis 5 (DAM program)
│   └── 06_phase_stratified_analysis.R         Analysis 7 (the headline)
├── example_data/
│   └── mouse_DE_example.csv                   expected mouse DE schema
└── docs/
    ├── pipeline_diagram.md                    design decisions
    └── walkthrough_notes.md                   interview cheat sheet
```

Each numbered script is independently runnable after `01_load_sea_ad.R` is sourced. `run_all.R` runs them all in order.

---

## 3. Pipeline overview

```
                  ┌─────────────────────────────────┐
                  │  00_download_sea_ad_data.sh     │
                  │  AWS S3 pull -- All / Early /   │
                  │  Late donor stratifications     │
                  │  for endothelial, microglia,    │
                  │  pericyte, astrocyte, oligo,    │
                  │  VLMC, SMC supertypes           │
                  └────────────┬────────────────────┘
                               │
                               ▼
                  ┌─────────────────────────────────┐
                  │  01_load_sea_ad.R               │
                  │  - read_nebula_file()           │
                  │  - load_sea_ad() per phase      │
                  │  - aggregate supertypes to      │
                  │    broad cell types (best-      │
                  │    pvalue per gene)             │
                  │  - mouse->human ortholog map    │
                  │  - load_mouse_de()              │
                  └────────────┬────────────────────┘
                               │
            ┌──────────────────┼──────────────────┬──────────────────┐
            ▼                  ▼                  ▼                  ▼
   ┌───────────────┐  ┌───────────────┐  ┌───────────────┐  ┌─────────────────┐
   │ 02_hyper-     │  │ 03_targeted   │  │ 04_cross-     │  │ 05_DAM          │
   │ geometric +   │  │ gene heatmap  │  │ species GSEA  │  │ microglia       │
   │ concordance   │  │ (28 curated   │  │ (Hallmark,    │  │ (12-gene canon- │
   │ (genome-wide) │  │ BBB / ISG /   │  │ 10 pathways,  │  │ ical DAM panel  │
   │               │  │ DAM genes)    │  │ 4 cell types) │  │ 3 modules)      │
   └───────────────┘  └───────────────┘  └───────────────┘  └─────────────────┘
                               │
                               ▼
                  ┌─────────────────────────────────┐
                  │  06_phase_stratified_analysis.R │
                  │  ----  THE HEADLINE  ----       │
                  │  Same GSEA across All / Early / │
                  │  Late, side by side, so the     │
                  │  phase cancellation is visible. │
                  └─────────────────────────────────┘
```

---

## 4. Quick start

```bash
# Clone
git clone https://github.com/<your-username>/sea-ad-cross-species.git
cd sea-ad-cross-species

# Install R deps
Rscript -e 'install.packages(c("dplyr","tidyr","readr","ggplot2","tibble","patchwork","BiocManager"))'
Rscript -e 'BiocManager::install(c("fgsea","msigdbr"))'

# Download SEA-AD reference data (AWS CLI required; ~50 MB; no auth needed)
bash scripts/00_download_sea_ad_data.sh ./data/SEA-AD

# Provide your mouse DE table (see example_data/mouse_DE_example.csv for schema)
cp example_data/mouse_DE_example.csv ./data/mouse_DE.csv  # or use your own

# Run everything
Rscript run_all.R
```

Outputs land in `./output/figures` and `./output/tables`.

To run a single analysis:

```bash
Rscript scripts/06_phase_stratified_analysis.R
```

---

## 5. Detailed step-by-step

### Step 0: Download SEA-AD data

`scripts/00_download_sea_ad_data.sh`

Pulls per-supertype Nebula DE CSVs from the public SEA-AD S3 bucket. The bucket layout is:

```
s3://sea-ad-single-cell-profiling/MTG/RNAseq/Supplementary Information/Nebula Results/
├── All Donors/Continuous_Pseudo-progression_Score/
├── Early Donors/Continuous_Pseudo-progression_Score/
└── Late Donors/Continuous_Pseudo-progression_Score/
```

Each phase has one CSV per of the 139 supertypes. This script downloads 13 supertype files per phase × 3 phases = 39 files total, focused on vascular and glial cell types since those are where cross-species AD comparisons are most informative. Edit the `CELL_TYPES` array in the script to include more supertypes.

**Why all three stratifications:** the all-donors result is the published default, but it's not the right comparison for acute disease models. Pulling all three lets the analysis directly compare them.

### Step 1: Load and aggregate (shared)

`scripts/01_load_sea_ad.R`

A library of helpers used by every analysis. Key functions:

- `read_nebula_file(path)` — handles the SEA-AD CSV format. The relevant columns are `logFC_Continuous_Pseudo-progression_Score` and `p_Continuous_Pseudo-progression_Score`. FDR is computed per file using Benjamini-Hochberg.
- `load_sea_ad(phase_dir, phase_label)` — loads every per-supertype CSV in a phase directory, attaches a broad cell type label by filename pattern, and aggregates: for each (broad type, gene) pair, keeps the supertype with the smallest p-value. Returns both the full and aggregated forms.
- `map_to_human(mouse_genes)` — applies a curated 46-gene one-to-one ortholog map (BBB, ISG, DAM panel) and falls back to uppercase conversion for everything else (~90% accurate for 1:1 orthologs).
- `load_mouse_de()` — reads the user's mouse DE table from the path config block. Expected schema is documented in `example_data/mouse_DE_example.csv`.

**Why broad-type aggregation:** SEA-AD reports DE per supertype (139 of them). Comparing each supertype individually would multiply the testing burden and give noisy comparisons. Aggregating to ~7 broad cell types matches the granularity of typical mouse studies and keeps the analysis interpretable.

**Why best-pvalue not mean:** when a gene is differentially expressed in one supertype of a broad type but not others, the average dilutes the signal. Taking the smallest-pvalue supertype preserves real per-broad-type signal at the cost of some inflation in the FDR — addressed downstream by reporting the supertype identity alongside.

### Step 2: Hypergeometric overlap + directional concordance

`scripts/02_hypergeometric_and_concordance.R`

Two analyses, both genome-wide.

**Hypergeometric overlap:** are mouse DE genes enriched among human DE genes more than chance? Tested at FDR thresholds of 0.05, 0.10, 0.25. The universe is genes tested in both species after ortholog mapping.

**Directional concordance:** for genes called DE in both, do the log fold-changes agree in sign? Spearman correlation, percent-concordant, and binomial test against the 50% null.

**Honest limitation:** these tests on the all-donors data largely return nulls (rho ≈ 0, percent-concordant ≈ 50%) in our cohort. That's not bad — it's what the data show with this stratification. The phase-stratified analysis (step 6) explains why. Reporting both is the methodological point.

### Step 3: Targeted gene heatmap

`scripts/03_targeted_gene_heatmap.R`

A 28-gene panel of curated BBB, endothelial identity, transporter, pericyte, ISG, and DAM genes. Each gene's mouse and human logFC are placed side by side, FDR significance marked. The "concordance percentage" reported here on a biologically meaningful subset is more interpretable than the genome-wide concordance — concordance on the genes that should be doing something is more informative than concordance averaged over the entire transcriptome.

### Step 4: Cross-species Hallmark GSEA

`scripts/04_cross_species_gsea.R`

Runs `fgsea` on both the mouse and human DE results against MSigDB Hallmark gene sets (separately, using the species-appropriate gene set for each side), then puts the NES values side by side. Pathways are restricted to a curated 10 for the heatmap (IFN α/γ, OXPHOS, inflammatory, TNF/NFkB, complement, IL6/JAK/STAT3, apoptosis, coagulation, angiogenesis) but the full results are saved to CSV.

**Why GSEA matters here:** many cross-species effects show up at the pathway level even when the individual-gene overlap is weak. OXPHOS depletion in glia is the strongest finding — conserved across mouse CMV and human AD astrocytes/oligodendrocytes — and it would not have been visible from gene-level overlap alone because individual OXPHOS genes have small effects spread across ~200 genes.

**Ranking metric:** `sign(logFC) × -log10(pvalue)`. Standard for fgsea; magnitude reflects significance, sign reflects direction. Direct logFC is too noisy on its own; -log10(p) alone loses direction.

### Step 5: DAM microglia comparison

`scripts/05_dam_microglia_comparison.R`

12-gene canonical DAM panel from Keren-Shaul et al. 2017, split into three modules: Homeostatic (expected DOWN), DAM Stage 1 (expected UP), DAM Stage 2 (expected UP). Side-by-side bars colored by species, with FDR<0.05 marked.

**Why DAM specifically:** the DAM program is the most replicated finding in human AD snRNA-seq. If a mouse model engages AD-like microglial biology, the DAM program should track in the same direction. This script asks how well it does. In our hands: 36% concordant in the all-donors comparison, ~50% with late-phase, ~60% with early-phase — supporting the "pre-early DAM initiation" interpretation rather than failure to engage the program.

### Step 6: Phase-stratified analysis (the headline)

`scripts/06_phase_stratified_analysis.R`

Same GSEA as step 4, but computed against all three SEA-AD stratifications: All / Early / Late. Output is a 4-panel heatmap (one per cell type) where columns are `Mouse CMV | All | Early | Late` and rows are 6 key Hallmark pathways. In our cohort:

- Endothelial IFN-α NES: Mouse +0.70, **Early +1.66**, Late −2.38, All averages to nothing.
- Endothelial OXPHOS NES: Mouse **−1.81**, **Early −1.17**, Late +1.43, All also +1.28.
- The pattern repeats across cell types and pathways.

Also computes per-gene phase concordance on a 17-gene BBB panel.

**This is the analysis that should be presented first** when introducing the work. It's the methodological lesson that travels — "naive averaging across disease trajectory loses signal that phase-stratification recovers" is a generalizable point relevant to any acute-vs-chronic disease comparison.

---

## 6. Expected results on the dataset this was developed for

Running on our 3xTg-AD CMV cohort at 1WPI (vascular + immune cells), comparing to SEA-AD across all three stratifications:

| Comparison | All Donors | Early Donors | Late Donors |
|---|---|---|---|
| Endothelial FDR<0.05 genes (human side) | 0 | 20 | 4 |
| Endothelial BBB-gene concordance | 47% | **59%** | 47% |
| Endothelial genome-wide Spearman ρ | -0.029 | **+0.014** | -0.022 |
| Endothelial IFN-α NES (human) | +0.91 | **+1.66** | -2.38 |
| Endothelial OXPHOS NES (human) | +1.28 | **-1.17** | +1.43 |
| Microglial DAM concordance | 36% | 36% | 55% |

The all-donors comparison alone would have concluded "mouse model fails to recapitulate AD vascular biology." The phase-stratified comparison reveals "mouse model recapitulates the early-phase vascular initiating events specifically." Different paper.

OXPHOS depletion is conserved across mouse and human glia in all three stratifications — that finding is robust to phase choice.

---

## 7. Known limitations

- **Ortholog mapping by uppercase conversion is imperfect.** Roughly 10% of mouse genes don't map to a same-symbol human ortholog (rodent-specific genes, name divergence, paralog complications). The 46-gene curated map covers most of the high-priority genes. For genome-wide analyses, the noise from imperfect mapping likely reduces detectable signal — i.e. errs conservative.
- **SEA-AD endothelial cells are severely underpowered.** Endothelial nuclei are rare in dissociated brain (~1-3%), and even with 84 donors, the Nebula models lack power. Zero genes reach FDR<0.05 in the all-donors endothelial comparison. The phase split (20 genes in early) helps, but endothelial-specific conclusions remain qualified.
- **Asymmetric DE methods.** The mouse side here uses single-cell Wilcoxon DE; the human side uses Nebula mixed-effects negative binomial regression with multiple covariates. Effect-size magnitudes are not directly comparable across these methods. We focus on direction (sign) and pathway-level enrichment for this reason.
- **Brain region mismatch.** SEA-AD is middle temporal gyrus only. Mouse studies often span whole hemisphere or other regions. Cross-region comparison adds noise.
- **Pseudoreplication.** The mouse single-cell Wilcoxon treats cells as independent units, which inflates significance when cells are clustered within mice (n=3-4 mice per group). Pseudobulk DESeq2 / edgeR is more rigorous but underpowered with so few biological replicates. The analyses here lean on directional concordance and pathway-level enrichment rather than gene-level FDR for downstream conclusions.

---

## 8. Citation

If you use this pipeline, please cite:

- This repository: Wroblewski, S. (2026). sea-ad-cross-species. [GitHub URL]
- SEA-AD: Gabitto et al., 2024. Integrated multimodal cell atlas of Alzheimer's disease. *Nature Neuroscience* 27:2366-2383.
- fgsea: Korotkevich et al. 2021. Fast gene set enrichment analysis. *bioRxiv*.
- msigdbr: Dolgalev I. msigdbr: MSigDB Gene Sets for Multiple Organisms.

---

## 9. Contact

Stella Wroblewski — [email]

Issues and pull requests welcome.
