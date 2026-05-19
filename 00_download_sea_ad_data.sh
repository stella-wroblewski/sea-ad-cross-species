#!/usr/bin/env bash
# =============================================================================
# 00_download_sea_ad_data.sh
# -----------------------------------------------------------------------------
# Downloads per-supertype Nebula differential expression results from the
# SEA-AD public AWS S3 bucket (Gabitto et al. 2024, Nature Neuroscience).
#
# Three "donor stratifications" are downloaded:
#   - All Donors (full CPS trajectory across 84 donors)
#   - Early Donors (early-phase pathology subset)
#   - Late Donors  (late-phase pathology subset)
#
# Only vascular and glial cell types are downloaded by default since this
# pipeline focuses on cross-species comparison with a CMV mouse model whose
# strongest signals are in those cell types. Edit the CELL_TYPES array to
# add additional supertypes.
#
# Requirements: AWS CLI v2.x. No AWS account needed — bucket is public
# (we use --no-sign-request).
#
# Usage: ./00_download_sea_ad_data.sh [output_directory]
# Default output: ./data/SEA-AD
#
# Author: Stella Wroblewski
# =============================================================================

set -e
set -u

OUT_DIR="${1:-./data/SEA-AD}"
mkdir -p "$OUT_DIR/All" "$OUT_DIR/Early" "$OUT_DIR/Late"

S3_ROOT="s3://sea-ad-single-cell-profiling/MTG/RNAseq/Supplementary Information/Nebula Results"

# Cell-type files to download. Each line is a single S3 filename relative to
# the donor-stratification subdir. Edit to include more supertypes if needed.
read -r -d '' CELL_TYPES <<'EOF' || true
Endo_Endo_1_across_Continuous_Pseudo-progression_Score_DE.csv
Endo_Endo_2_across_Continuous_Pseudo-progression_Score_DE.csv
Endo_Endo_3_across_Continuous_Pseudo-progression_Score_DE.csv
Micro-PVM_Micro-PVM_1_across_Continuous_Pseudo-progression_Score_DE.csv
Micro-PVM_Micro-PVM_2_across_Continuous_Pseudo-progression_Score_DE.csv
VLMC_Pericyte_1_across_Continuous_Pseudo-progression_Score_DE.csv
VLMC_Pericyte_2-SEAAD_across_Continuous_Pseudo-progression_Score_DE.csv
VLMC_VLMC_1_across_Continuous_Pseudo-progression_Score_DE.csv
VLMC_SMC-SEAAD_across_Continuous_Pseudo-progression_Score_DE.csv
Astro_Astro_1_across_Continuous_Pseudo-progression_Score_DE.csv
Astro_Astro_2_across_Continuous_Pseudo-progression_Score_DE.csv
Oligo_Oligo_1_across_Continuous_Pseudo-progression_Score_DE.csv
Oligo_Oligo_2_across_Continuous_Pseudo-progression_Score_DE.csv
EOF

# Check AWS CLI
if ! command -v aws &>/dev/null; then
  echo "ERROR: AWS CLI not found. Install from https://aws.amazon.com/cli/"
  exit 1
fi

echo "Downloading SEA-AD Nebula DE results to: $OUT_DIR"
echo "Source: $S3_ROOT"
echo ""

download_phase() {
  local phase_label="$1"   # All | Early | Late
  local phase_subdir="$2"  # "All Donors" | "Early Donors" | "Late Donors"
  local dest="$OUT_DIR/$phase_label"

  echo "=== $phase_label Donors ==="
  while IFS= read -r f; do
    [ -z "$f" ] && continue
    if [ -f "$dest/$f" ]; then
      echo "  [skip] $f already present"
      continue
    fi
    local src="$S3_ROOT/$phase_subdir/Continuous_Pseudo-progression_Score/$f"
    echo "  downloading $f"
    aws s3 cp --no-sign-request "$src" "$dest/" --quiet || \
      echo "  [WARN] failed to fetch $f (may not exist for this phase)"
  done <<< "$CELL_TYPES"
  echo ""
}

download_phase "All"   "All Donors"
download_phase "Early" "Early Donors"
download_phase "Late"  "Late Donors"

echo "Done."
echo "Files in $OUT_DIR:"
find "$OUT_DIR" -name "*.csv" | sort
