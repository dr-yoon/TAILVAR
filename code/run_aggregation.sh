#!/bin/bash
set -euo pipefail  # Exit on error, treat unset variables as an error, and pipe failures propagate

tango="/path_to/tango_2.3.1/tango_x86_64_release" # Path to TANGO
WD=$1
file_name=$2
SEQDIR="${WD}/SEQ_${file_name}"
seq_input="${WD}/${file_name}_seq_input.txt"
tango_result="${WD}/${file_name}_TANGO_score.tsv"
canya_result="${WD}/${file_name}_CANYA_score.tsv"

mkdir -p ${SEQDIR}
cd ${SEQDIR} || exit 1

# run TANGO for each peptide
while IFS=$'\t' read -r Name Sequence; do
    if [[ -z "$Sequence" ]]; then
        echo "Empty sequence for ${Name}, skipping."
        continue
    fi
    echo "Running TANGO for $Name"
    "$tango" "$Name" ct="N" nt="N" ph="7.4" te="310" io="0.15" seq="$(echo "$Sequence" | tr -d '\r\n')"
done < "$seq_input"

echo -e "Variant\tTANGO" > "$tango_result" 

# combine TANGO scores into one file
for file in "${SEQDIR}"/*.txt; do
  [[ -e "$file" ]] || continue
        
  variant=$(basename "$file" .txt)
  max_score=$(awk -F'\t' 'NR > 1 { gsub(/^[ \t]+/, "", $6); if ($6 > max) max = $6 } END { print (max == "") ? "0.000" : max }' "$file")

  printf '%s\t%s\n' "$variant" "$max_score" >> "$tango_result"
done

# run CANYA
conda activate canyaenv
canya --input ${seq_input} --output ${canya_result}
