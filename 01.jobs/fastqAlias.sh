
#!/bin/bash

# Make sure you're in the directory with the Ca*.fastq.gz files
cd ~/workspace/rnaseq/data/Pi/Ca || exit 1

# Remove previous symlinks if re-running
##rm -f Ca_1.*.fastq.gz Ca_2.*.fastq.gz

# Initialize counter
i=1

# Loop through R1 files only, in sorted order
for r1 in $(ls Ca*_R1.fastq.gz | sort); do
    # Derive R2 filename
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"

    # Check both files exist
    if [[ -f "$r1" && -f "$r2" ]]; then
        ln -s "$r1" "Ca_1.${i}.fastq.gz"
        ln -s "$r2" "Ca_2.${i}.fastq.gz"
        echo "Linked $r1 → Ca_1.${i}.fastq.gz"
        echo "Linked $r2 → Ca_2.${i}.fastq.gz"
        ((i++))
    else
        echo "⚠️ Skipping $r1 or $r2 — missing file."
    fi
done

