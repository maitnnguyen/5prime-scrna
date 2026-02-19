#!/usr/bin/env bash
# ============================================================
#  bin/fixed_bidir_enhancers.sh
#  Bidirectional enhancer detection script
#  Adapted from FANTOM5 consortium (Andersson et al., Nature 2014)
#  and ReapTEC (Oguchi et al., Science 2024)
#  Original: github.com/anderssonrobin/enhancers
#  ReapTEC version: github.com/MurakawaLab/ReapTEC
# ============================================================
#
#  USAGE:
#    fixed_bidir_enhancers.sh <sense.bed> <antisense.bed> <gap_bp> <out_prefix>
#
#  INPUTS:
#    sense.bed      - BED file of sense strand TSS peaks (+)
#    antisense.bed  - BED file of antisense strand TSS peaks (-)
#    gap_bp         - Max gap between sense/antisense pair (default: 500)
#    out_prefix     - Output file prefix
#
#  OUTPUT:
#    <out_prefix>.bidirectional_enhancers.bed
#      Format: chr  start  end  name  score  strand
#              (enhancer center ± window, bidirectional score)
#
# ============================================================

set -euo pipefail

SENSE_BED=${1:?"ERROR: sense BED required as arg 1"}
ANTI_BED=${2:?"ERROR: antisense BED required as arg 2"}
GAP=${3:-500}
PREFIX=${4:-"bidirectional"}
WINDOW=${5:-200}
MIN_CPM=${6:-2}

echo "[fixed_bidir_enhancers] ========================="
echo "[fixed_bidir_enhancers] Sense peaks:     $SENSE_BED"
echo "[fixed_bidir_enhancers] Antisense peaks: $ANTI_BED"
echo "[fixed_bidir_enhancers] Max gap:         ${GAP} bp"
echo "[fixed_bidir_enhancers] Window:          ±${WINDOW} bp"
echo "[fixed_bidir_enhancers] Min log2CPM:     ${MIN_CPM}"
echo "[fixed_bidir_enhancers] Output prefix:   $PREFIX"
echo "[fixed_bidir_enhancers] ========================="

# ── STEP 1: Validate inputs ──────────────────────────────────
if [ ! -f "$SENSE_BED" ]; then
    echo "ERROR: Sense BED file not found: $SENSE_BED" >&2; exit 1
fi
if [ ! -f "$ANTI_BED" ]; then
    echo "ERROR: Antisense BED file not found: $ANTI_BED" >&2; exit 1
fi

SENSE_N=$(wc -l < "$SENSE_BED")
ANTI_N=$(wc -l < "$ANTI_BED")
echo "[fixed_bidir_enhancers] Input: $SENSE_N sense peaks | $ANTI_N antisense peaks"

# ── STEP 2: Merge CTSS within 10bp on same strand ────────────
# This clusters individual TSS positions into broader peaks
# 10bp merge window is the ReapTEC default from fixed_bidir_enhancers_10bp.sh

echo "[fixed_bidir_enhancers] Step 2: Merging CTSS within 10bp..."

bedtools merge -d 10 -i "$SENSE_BED" \
    -c 5 -o sum \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"sense_"NR,$4,"+"}' \
    | sort -k1,1 -k2,2n > tmp_sense_merged.bed

bedtools merge -d 10 -i "$ANTI_BED" \
    -c 5 -o sum \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"anti_"NR,$4,"-"}' \
    | sort -k1,1 -k2,2n > tmp_anti_merged.bed

SENSE_MERGED=$(wc -l < tmp_sense_merged.bed)
ANTI_MERGED=$(wc -l < tmp_anti_merged.bed)
echo "[fixed_bidir_enhancers] After 10bp merge: $SENSE_MERGED sense | $ANTI_MERGED antisense"

# ── STEP 3: Find bidirectional pairs within gap ───────────────
# For each sense peak, find antisense peaks within GAP bp
# The pair must be on opposite strands (enforced by -Sm flag)

echo "[fixed_bidir_enhancers] Step 3: Pairing sense/antisense peaks within ${GAP}bp..."

bedtools window \
    -a tmp_sense_merged.bed \
    -b tmp_anti_merged.bed \
    -w "$GAP" \
    -Sm \
    | awk -v win="$WINDOW" 'BEGIN{OFS="\t"}{
        # Sense peak center
        s_center = int(($2 + $3) / 2);
        # Antisense peak center
        a_center = int(($8 + $9) / 2);
        # Enhancer center = midpoint
        enh_center = int((s_center + a_center) / 2);
        enh_start  = enh_center - win;
        enh_end    = enh_center + win;
        if (enh_start < 0) enh_start = 0;
        # Bidirectional score = sum of both strand signals
        bidir_score = $5 + $11;
        distance    = a_center - s_center;
        print $1, enh_start, enh_end, \
              "btcEnh_"NR, bidir_score, ".", \
              s_center, a_center, $5, $11, distance
    }' \
    | sort -k1,1 -k2,2n > tmp_bidir_raw.bed

RAW_N=$(wc -l < tmp_bidir_raw.bed)
echo "[fixed_bidir_enhancers] Raw bidirectional pairs: $RAW_N"

# ── STEP 4: Merge overlapping enhancer candidates ────────────
echo "[fixed_bidir_enhancers] Step 4: Merging overlapping candidates..."

cut -f1-6 tmp_bidir_raw.bed \
    | bedtools merge -i stdin -c 4,5 -o first,sum \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"btcEnh_"NR,$5,"."}' \
    | sort -k1,1 -k2,2n > tmp_merged_enhancers.bed

MERGED_N=$(wc -l < tmp_merged_enhancers.bed)
echo "[fixed_bidir_enhancers] After merge: $MERGED_N candidates"

# ── STEP 5: Apply CPM score filter ───────────────────────────
echo "[fixed_bidir_enhancers] Step 5: Applying log2CPM >= ${MIN_CPM} filter..."

awk -v min_cpm="$MIN_CPM" '
    $5 >= min_cpm { print }
' tmp_merged_enhancers.bed > "${PREFIX}.bidirectional_enhancers.bed"

FINAL_N=$(wc -l < "${PREFIX}.bidirectional_enhancers.bed")
echo "[fixed_bidir_enhancers] Final enhancers (log2CPM >= ${MIN_CPM}): $FINAL_N"

# ── STEP 6: Summary statistics ───────────────────────────────
echo ""
echo "[fixed_bidir_enhancers] ====== SUMMARY ======"
echo "[fixed_bidir_enhancers] Sense peaks input:         $SENSE_N"
echo "[fixed_bidir_enhancers] Antisense peaks input:     $ANTI_N"
echo "[fixed_bidir_enhancers] After 10bp merge (sense):  $SENSE_MERGED"
echo "[fixed_bidir_enhancers] After 10bp merge (anti):   $ANTI_MERGED"
echo "[fixed_bidir_enhancers] Raw bidirectional pairs:   $RAW_N"
echo "[fixed_bidir_enhancers] After candidate merge:     $MERGED_N"
echo "[fixed_bidir_enhancers] Final btcEnhs:             $FINAL_N"
echo "[fixed_bidir_enhancers] Output: ${PREFIX}.bidirectional_enhancers.bed"
echo "[fixed_bidir_enhancers] ====================="

# ── Cleanup temp files ────────────────────────────────────────
rm -f tmp_sense_merged.bed tmp_anti_merged.bed tmp_bidir_raw.bed tmp_merged_enhancers.bed

echo "[fixed_bidir_enhancers] Done."
