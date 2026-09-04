#!/usr/bin/env bash

# Driver supplied for the Apennine brown bear GONE analysis.
# Usage: bash script_GONE.sh FILE WORKING_DIRECTORY

FILE=$1
WD=$2

source INPUT_PARAMETERS_FILE

if [ -f "OUTPUT_$FILE" ]; then rm "OUTPUT_$FILE"; fi
if [ -f "Ne_$FILE" ]; then rm "Ne_$FILE"; fi
if [ -d "TEMPORARY_FILES" ]; then rm -r "TEMPORARY_FILES"; fi
mkdir TEMPORARY_FILES

cp "$WD/$FILE.map" data.map
cp "$WD/$FILE.ped" data.ped

tr '\t' ' ' < data.map > KK1
cut -d ' ' -f1 < KK1 > KK2
grep -w "" -c data.ped > NIND
tail -n 1 data.map | tr '\t' ' ' | cut -d ' ' -f1 > NCHR

SAM=$(grep -w "" -c data.ped)
NCHR=$(tail -n 1 data.map | tr '\t' ' ' | cut -d ' ' -f1)

for ((i=1; i<=NCHR; i++)); do
  grep -wc "$i" < KK2 > "NCHR$i"
done

if [ -f "SNP_CHROM" ]; then rm SNP_CHROM; fi
for ((i=1; i<=NCHR; i++)); do cat "NCHR$i" >> SNP_CHROM; done
rm KK1 KK2

echo "DIVIDE .ped AND .map FILES IN CHROMOSOMES"
echo "DIVIDE .ped AND .map FILES IN CHROMOSOMES" > timefile
num=$RANDOM
echo "$num" > seedfile

./PROGRAMMES/MANAGE_CHROMOSOMES2 >> out <<EOF
-99
$maxNSNP
EOF

rm NCHR*
rm NIND
rm SNP_CHROM

if [ "$maxNCHROM" != -99 ]; then NCHR=$maxNCHROM; fi

echo "RUNNING ANALYSIS OF CHROMOSOMES ..."
echo "RUNNING ANALYSIS OF CHROMOSOMES" >> timefile

options_for_LD="$SAM $MAF $PHASE $NGEN $NBIN $ZERO $DIST $cMMb"
if [ "$threads" -eq -99 ]; then threads=$(getconf _NPROCESSORS_ONLN); fi

START=$(date +%s)
cp chromosome* TEMPORARY_FILES/

for ((n=1; n<=NCHR; n++)); do echo "$n"; done | \
  xargs -I % -P "$threads" bash -c "./PROGRAMMES/LD_SNP_REAL3 % $options_for_LD"

END=$(date +%s)
DIFF=$((END - START))
echo "CHROMOSOME ANALYSES took $DIFF seconds"
echo "CHROMOSOME ANALYSES took $DIFF seconds" >> timefile

for ((n=1; n<=NCHR; n++)); do
  cat "outfileLD$n" >> CHROM
  echo "CHROMOSOME $n" >> OUTPUT
  sed '2,3d' "outfileLD$n" > temp
  mv temp "outfileLD$n"
  cat "parameters$n" >> OUTPUT
done

mv outfileLD* TEMPORARY_FILES/
rm parameters*

./PROGRAMMES/SUMM_REP_CHROM3 >> out <<EOF
$NGEN   NGEN
$NBIN   NBIN
$NCHR   NCHR
EOF

mv chrom* TEMPORARY_FILES/

echo "TOTAL NUMBER OF SNPs" >> "OUTPUT_$FILE"
cat nsnp >> "OUTPUT_$FILE"
echo -e "\n" >> "OUTPUT_$FILE"
echo "HARDY-WEINBERG DEVIATION" >> "OUTPUT_$FILE"
cat outfileHWD >> "OUTPUT_$FILE"
echo -e "\n" >> "OUTPUT_$FILE"
cat OUTPUT >> "OUTPUT_$FILE"
echo -e "\n" >> "OUTPUT_$FILE"
echo "INPUT FOR GONE" >> "OUTPUT_$FILE"
echo -e "\n" >> "OUTPUT_$FILE"
cat outfileLD >> "OUTPUT_$FILE"

rm nsnp OUTPUT CHROM

echo "Running GONE"
echo "Running GONE" >> timefile
START=$(date +%s)
./PROGRAMMES/GONEparallel.sh -hc "$hc" outfileLD "$REPS"
END=$(date +%s)
DIFF=$((END - START))
echo "GONE run took $DIFF seconds"
echo "GONE run took $DIFF seconds" >> timefile
echo "END OF ANALYSES"
echo "END OF ANALYSES" >> timefile

mv outfileLD_Ne_estimates "Output_Ne_$FILE"
mv outfileLD_d2_sample "Output_d2_$FILE"
rm outfileLD data.ped data.map out
mv outfileLD_TEMP TEMPORARY_FILES/

mv *"$FILE"* "$WD"
mv timefile seedfile TEMPORARY_FILES/ outfileHWD "$WD"
