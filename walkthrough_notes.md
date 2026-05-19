# Code walkthrough notes (interview prep) -- SEA-AD repo

This is a backup option for the code walkthrough, or to mention alongside the Tangram repo if Mike asks about Alzheimer's-specific applications. The narrative is different from the Tangram walkthrough: that one is about ML methods, this one is about a methodological insight.

## The story this code tells

The headline is the **phase-cancellation insight**. Most people who download SEA-AD use the all-donors Nebula coefficients because they're the published default and the largest sample size. For chronic disease that's fine. For acute disease models — like our peripheral CMV infection at one week — the all-donors average dilutes early-phase signals through cancellation with late-phase reversals. Acute models should be compared to the early-phase stratification, not the all-donors aggregate.

I figured this out the hard way: ran the obvious all-donors comparison, got nothing significant, almost concluded the mouse model didn't engage AD biology. Re-read the Gabitto paper, noticed the two-phase trajectory, re-downloaded the phase-stratified files, saw a huge signal in early-phase that was hidden by cancellation in all-donors.

## Suggested walkthrough order (~20 min)

### 1. Open with the problem in the README (3 min)
- Section 1: "why this pipeline exists"
- The table showing early vs late phase characteristics
- Quick story: I almost wrote the wrong conclusion

### 2. The shared loader: `01_load_sea_ad.R` (3 min)
- `read_nebula_file()` — points at the relevant SEA-AD CSV columns
- `load_sea_ad()` — broad-cell-type aggregation rule
- The ortholog mapping strategy (curated for critical genes, uppercase fallback for everything else)
- This is reused everywhere; one place to maintain

### 3. The naive analyses (3 min)
- Quick walkthrough of `02_hypergeometric_and_concordance.R`
- Mention these largely return null on all-donors
- Important: I report the nulls honestly. The trap is presenting only the analyses that "worked"

### 4. The pathway-level approach (3 min)
- `04_cross_species_gsea.R`
- Why GSEA picks up effects that gene-level overlap misses
- Pathway integration averages over the ortholog noise
- OXPHOS depletion in glia is the conserved finding

### 5. The headline: phase-stratified (5 min)
- `06_phase_stratified_analysis.R` — the script to spend time on
- Walk through the heatmap structure: Mouse | All | Early | Late × cell types
- Point at the IFN-α row: Mouse +0.70, Early +1.66, Late −2.38, All +0.91
- The "All" column averages early and late and gets nothing
- The "Early" column reveals the real cross-species concordance

### 6. The DAM caveat (2 min)
- `05_dam_microglia_comparison.R`
- Honest finding: DAM concordance is weaker (36-55%)
- Biologically sensible interpretation: 1WPI is "pre-early" — microglia are just starting the transition while early-phase AD microglia have already completed it
- Important to flag what doesn't work, not just what does

### 7. Limitations (1 min)
- README section 7. Be honest:
  - Endothelial cells underpowered in SEA-AD
  - Pseudoreplication concern in mouse Wilcoxon DE
  - Brain region mismatch (MTG vs whole hemisphere)
  - Ortholog mapping imperfect

## Anticipated questions

**"Why didn't you use Supplementary Table 5 from the Gabitto paper instead of the AWS bucket?"**
- Supp Table 5 is all-donors only. The phase stratifications are on the S3 bucket. So either way, the bucket is the actual data source for this analysis.

**"Why FDR per file rather than across all files?"**
- This matches what SEA-AD does internally and preserves comparability with their published results. A global FDR across all supertypes would be more conservative but would also conflate cell-type-specific significance with overall significance.

**"Best-pvalue aggregation across supertypes -- doesn't that inflate the FDR?"**
- Yes, somewhat. The alternatives (mean logFC, max effect size, Stouffer's meta-analysis) all have trade-offs. The best-pvalue rule preserves real per-broad-type biology at the cost of mild p-value inflation. Reporting which supertype contributed the smallest p-value alongside lets readers calibrate.

**"You're treating mouse Wilcoxon and human Nebula effect sizes as comparable. They're not."**
- Right. The text and figures emphasize direction and pathway-level enrichment for this reason. The hypergeometric tests don't require magnitude comparability — only thresholded membership. The GSEA uses ranks. The targeted heatmap shows magnitudes but they're meant to be read within-species.

**"What about ADNC stratification rather than the binary early/late split?"**
- SEA-AD also provides DE results stratified by ADNC (no-AD / low / intermediate / high). I haven't done that comparison yet — it would give finer resolution and likely show the same gradient. Listed in future directions.

**"Why R for this and not Python?"**
- The SEA-AD release is R-friendly, `msigdbr` is the cleanest MSigDB interface available, and `fgsea` is the fastest GSEA implementation. There's no reason to move to Python except interop with downstream ML, which isn't part of this pipeline.

**"How would you extend this to other mouse models?"**
- The mouse DE file is the only model-specific input. Swap it for 5xFAD, APP/PS1, Tg2576, or any other model's DE table and re-run. The phase-stratified comparison would let us ask whether different mouse models recapitulate early vs late AD differently — a generalization that would itself be publishable.

**"What was the AI / ML angle here?"** (if Mike asks since this is for an AI/Alzheimer's role)
- This specific repo is more bioinformatics than ML, but it sits inside a three-layer framework with the Kosmos AI platform (Future House / Edison Scientific) as the discovery layer. Kosmos ran ~200 AI agents over the dataset and converged on vascular cells as the responders; the manual validation pipeline and this cross-species comparison are the verification layers. The AI angle is the validation framework for AI-generated hypotheses, not the cross-species code itself.

## Honesty notes — things you should NOT oversell

- **The data is from SEA-AD; I didn't generate it.** My contribution is the cross-species comparison pipeline and the phase-stratification insight.
- **The phase-cancellation point isn't new in the SEA-AD paper.** Gabitto et al. document the two phases. What I did was operationalize it as a comparison strategy for an acute mouse model.
- **The endothelial findings are qualified by SEA-AD's underpowered endothelial data.** Even the early-phase 20-gene FDR<0.05 result is small in absolute terms.
- **The microglial DAM comparison doesn't work as well as the vascular comparison.** That's in the analysis. Don't bury it.
