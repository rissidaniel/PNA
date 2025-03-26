
!/bin/bash

echo "generating_primer_candidates .... no trimming"

# Loop through -d 1 and 2, and -l 12 to 24 for untrimmed alignments
for d in 1 2; do
    for l in {12..24}; do
        perl /data/daniel/pr2/DEGEPRIME/DegePrime.pl \\
            -i clustalo_targerts.pcr.fasta_rep_set.fasta \\
            -d $d -l $l \\
            -o output_file_${d}_${l}.txt
    done
done

# Run trimming at different thresholds and repeat primer generation
for trim in 70 80 90; do
    perl /data/daniel/pr2/DEGEPRIME/TrimAlignment.pl \\
        -i clustalo_targerts.pcr.fasta_rep_set.fasta \\
        -min 0.$trim \\
        -o trimmed_align_file_$trim

    for d in 1 2; do
        for l in {12..24}; do
            perl /data/daniel/pr2/DEGEPRIME/DegePrime.pl \\
                -i trimmed_align_file_$trim \\
                -d $d -l $l \\
                -o output_file_${trim}_${d}_${l}.txt
        done
    done
done

# Combine and extract primer sequences
cat output_file_*.txt > tmp_primer_candidates.txt
awk '{print $7}' tmp_primer_candidates.txt | sort | uniq -u > primer_candidates.txt

# Format FASTA
sed -i 's/PrimerSe//g' primer_candidates.txt
sed -i '/^$/d' primer_candidates.txt
sed -i 's/^/>/' primer_candidates.txt

# Split into individual files
awk '{ if (substr($0,1,1)==">") {filename=(substr($0,2) ".text")} print $0 > filename }' primer_candidates.txt
sed -i 's/>/forward     /g' *.text

# Mothur analysis
for f in *.text; do
    SAMPLE=$(basename "${f%%.text}")
    mothur "#pcr.seqs(fasta=targets.pcr.fasta, oligos=${SAMPLE}.text, processors=40);"
    mv targets.pcr.pcr.fasta targets.pcr.${SAMPLE}.fasta
done

# Summarize results
grep -c "Fungi" targets.pcr.*.fasta > results_fungi
grep -c ">" targets.pcr.*.fasta > results_all

sed -i '1 i\\ncandidate:fungi' results_fungi
sed -i '1 i\\ncandidate:total' results_all

sed -i 's/:/\\t/' results_all
sed -i 's/:/\\t/' results_fungi

paste results_all results_fungi | awk '{print $1,$2,$4}' > results_comparative

echo "results in file: results_comparative"
echo "the end!!!"
