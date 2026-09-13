#!/usr/bin/env bash
# Instructor sensitivity analysis: 36 scaffolds x 2 populations x 2 -u values.
# Run from the Day 1 course directory, after activating bear-ne-practical.
set -euo pipefail

SOURCE_VCF=/jarvis/scratch/usr/biello/test/day01/data/Ursus.allchr.hq.snp.masked.filt.GQ10.noallmiss.norepeat.highmap.noHetExcess.vcf.gz
GONE2_BIN="$PWD/software/GONE2/gone2"
ROOT="$PWD/results/gone2_per_scaffold"
INPUT="$ROOT/input"
OUTPUT="$ROOT/output"
LOGS="$ROOT/logs"
MANIFEST="$ROOT/run_status.tsv"
mkdir -p "$INPUT" "$OUTPUT" "$LOGS"

for executable in bcftools plink; do command -v "$executable" >/dev/null || { echo "Missing: $executable" >&2; exit 1; }; done
[[ -x "$GONE2_BIN" ]] || { echo "Missing GONE2: $GONE2_BIN" >&2; exit 1; }
[[ -s "$SOURCE_VCF" ]] || { echo "Missing source VCF: $SOURCE_VCF" >&2; exit 1; }
for sample_file in data/apennine.samples data/slovak.samples; do
  [[ -s "$sample_file" ]] || { echo "Missing sample list: $sample_file" >&2; exit 1; }
done
printf 'Scaffold\tPopulation\tu\tMarkers\tSpan_bp\tStatus\n' > "$MANIFEST"

for number in $(seq 1 36); do
  chrom="Scaffold_${number}"
  short="s${number}"
  common="$INPUT/UrArMa_18i_${short}.vcf.gz"
  echo "Preparing $chrom"

  # Apply the site filters once, before splitting populations.
  if ! bcftools view -r "$chrom" \
      -S <(awk '{print $1}' data/apennine.samples data/slovak.samples | sort -u) \
      -Ou "$SOURCE_VCF" |
      bcftools view -v snps -m2 -M2 -e 'GT="mis"' \
        -Oz -o "$common"; then
    echo "Cannot extract $chrom" >&2
    exit 1
  fi
  bcftools index -f "$common"

  for population in ABB SBB; do
    if [[ "$population" == ABB ]]; then samples=data/apennine.samples; else samples=data/slovak.samples; fi
    pop_vcf="$INPUT/${population}_${short}.vcf.gz"
    prefix="$INPUT/${population}_${short}"
    temporary="${prefix}_temporary"

    bcftools view -S "$samples" -Oz -o "$pop_vcf" "$common"
    bcftools index -f "$pop_vcf"
    plink --vcf "$pop_vcf" --double-id --allow-extra-chr \
      --snps-only just-acgt --biallelic-only strict \
      --recode --out "$temporary" > "$LOGS/${population}_${short}_plink.log" 2>&1
    cp "${temporary}.ped" "${prefix}.ped"
    # GONE2 sees a single chromosome numbered 1; preserve the physical positions.
    awk 'BEGIN {OFS="\t"} {$1=1; print $1,$2,$3,$4}' \
      "${temporary}.map" > "${prefix}.map"

    markers=$(wc -l < "${prefix}.map" | tr -d ' ')
    span=$(awk 'NR==1{min=$4; max=$4} {if ($4<min) min=$4; if ($4>max) max=$4} END{if (NR) print max-min; else print 0}' "${prefix}.map")
    for u in 0.02 0.05; do
      tag=${u/./}  # 0.02 -> 002; 0.05 -> 005
      out="$OUTPUT/${population}_${short}_u${tag}"
      log="$LOGS/${population}_${short}_u${tag}.log"
      if (( markers < 10 || span <= 20000000 )); then
        printf '%s\t%s\t%s\t%s\t%s\tSKIP_short_or_sparse\n' "$chrom" "$population" "$u" "$markers" "$span" >> "$MANIFEST"
        continue
      fi
      if "$GONE2_BIN" -g 0 -u "$u" -r 1 -t 2 -S 1 -E \
          -o "$out" "${prefix}.ped" > "$log" 2>&1 && [[ -s "${out}_GONE2_Ne" ]]; then
        status=OK
      else
        status=FAILED
        echo "GONE2 failed: $chrom $population -u $u (see $log)" >&2
      fi
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$chrom" "$population" "$u" "$markers" "$span" "$status" >> "$MANIFEST"
    done
  done
done

echo "Completed. Check $MANIFEST before plotting."
