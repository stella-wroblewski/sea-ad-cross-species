# Pipeline architecture and design decisions

This document records the technical decisions behind the SEA-AD comparison pipeline.

## High-level data flow

```
                        ┌────────────────────────────┐
                        │  Mouse single-cell DE      │
                        │  (Wilcoxon, per cell class │
                        │   per timepoint)           │
                        └─────────────┬──────────────┘
                                      │ (one CSV)
                                      ▼
                        ┌────────────────────────────┐
                        │  load_mouse_de()           │
                        │  + map_to_human() ortholog │
                        └─────────────┬──────────────┘
                                      │
                                      ▼
                        ┌────────────────────────────┐
                        │  Joined per-gene table:    │
                        │  (mouse, human, logFC×2,   │
                        │   pvalue×2, FDR×2)         │
                        └─────────────┬──────────────┘
                                      │
            ┌─────────────────────────┼─────────────────────────┐
            ▼                         ▼                         ▼
   ┌────────────────┐        ┌────────────────┐        ┌────────────────┐
   │ phyper test    │        │ fgsea Hallmark │        │ Targeted gene  │
   │ (overlap)      │        │ (10 pathways)  │        │ panels (BBB,   │
   │                │        │                │        │ DAM)           │
   └────────────────┘        └────────────────┘        └────────────────┘

  All of the above run for All, Early, AND Late stratifications.
```

## Why R for everything

The SEA-AD release is R-friendly: published Nebula DE results are CSVs, MSigDB has a high-quality R interface via `msigdbr`, and `fgsea` is the fastest published GSEA implementation. Going to Python would add an interop layer with no real benefit.

## Key choice: aggregate supertypes to broad cell types

SEA-AD reports DE per supertype (139 of them, including 3 endothelial subtypes, 2 microglia, etc.). Two reasons we aggregate to broad cell types instead of comparing per-supertype:

1. **Mouse studies typically don't have 139-supertype resolution.** Our mouse cell classes are Allen taxonomy level-1 broad classes. The comparison only works at the resolution that exists on both sides.
2. **Per-supertype comparison would multiply the testing burden.** 139 × 4 mouse cell classes = 556 tests with no signal-aware multiple-testing control. Aggregating first is cleaner.

The aggregation rule is "smallest p-value supertype per gene." Alternatives:
- **Mean logFC:** dilutes signal when DE is supertype-specific.
- **Sum:** sensitive to outlier supertypes.
- **Meta-analysis (Stouffer):** more principled but requires per-supertype N which isn't always reported cleanly. Worth a future ablation.

## Key choice: phase-stratified comparison

The published default is to use the all-donors CPS coefficient. We instead compute three stratifications side-by-side because of the phase-cancellation problem described in the README.

Trade-off: smaller sample sizes per phase reduce statistical power within each stratification. We accept this because the alternative (averaging) destroys real biological signal in an opaque way. Reporting all three lets the reader see when this matters.

## Key choice: ortholog mapping by curated map + uppercase fallback

Two-tier approach:
1. **Curated 46-gene map** for the genes most relevant to the BBB / ISG / DAM comparison. These are the genes I actually need to be right.
2. **Uppercase conversion** for everything else. About 90% accurate for 1:1 mouse-human orthologs by gene symbol.

Alternatives:
- **HomoloGene / Ensembl Compara:** more thorough but adds a dependency and runtime hit. Worth doing for a publication-grade analysis; overkill for an exploratory pipeline.
- **biomaRt:** API-dependent, often slow or down.

The uppercase fallback errs conservative — wrong mappings just add noise to the genome-wide concordance, which biases toward null rather than toward false positive.

## Key choice: pathway-level GSEA over gene-level analysis as the primary readout

Gene-level concordance is noisy across species. Pathway-level GSEA integrates information across hundreds of genes per pathway and is less sensitive to:
- Ortholog mapping errors
- Per-gene technical noise
- The mouse-Wilcoxon vs human-Nebula effect-size scale mismatch

The strongest finding in our analysis (conserved OXPHOS depletion in glia) emerged from GSEA, not from gene-level overlap.

## Key choice: matched cell-type pairs are hard-coded

The `CELL_PAIRS` list is:
```r
list(
  list(mouse = "Vascular",   human = "Endothelial"),
  list(mouse = "Immune",     human = "Microglia"),
  list(mouse = "Astro-Epen", human = "Astrocyte"),
  list(mouse = "OPC-Oligo",  human = "Oligodendrocyte")
)
```

These are the obvious functional pairs. Pericytes show up in SEA-AD as a separate broad type but our mouse "Vascular" class lumps endothelial + pericyte cells together, so we don't split it. If someone has finer-resolution mouse data, edit the list.

## Things this pipeline does NOT do

Kept simple to be a transparent baseline.

- **No pseudobulk DE on the mouse side.** Wilcoxon single-cell DE is used. We've separately discussed pseudobulk (DESeq2 / edgeR) as more rigorous, but with n=3-4 mice per group, pseudobulk methods have insufficient residual degrees of freedom to estimate dispersion. The pragmatic choice is single-cell Wilcoxon with the caveat that pseudoreplication inflates significance.
- **No formal cross-species effect-size harmonization.** The mouse logFC and the SEA-AD CPS coefficient are not on the same scale. We focus on sign (direction) and pathway enrichment, which are scale-invariant.
- **No per-supertype comparisons in the main analyses.** Could be added; complicates interpretation.
- **No incorporation of GWAS evidence or AD risk genes.** A natural extension — overlap our DE with AD risk loci. Worth adding.

## Future directions

- Add ADNC-stratified analysis (low/intermediate/high) as a finer phase resolution.
- Pseudobulk DE on the mouse side with appropriate covariates.
- Cross-comparison with other mouse AD models (5xFAD, APP/PS1, Tg2576) to see whether the early-phase concordance is CMV-specific or generalizable.
- Use scib-style metrics on the aligned cross-species datasets if the goal is a quantitative cross-species integration rather than a comparison.
