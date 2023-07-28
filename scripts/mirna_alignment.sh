#!/usr/bin/env bash

# Format miRBase hairpin file
if [ -f hairpin.fa ]; then
    :
else
    wget https://www.mirbase.org/download_file/hairpin.fa
fi
sed '#^[^>]#s#[^AUGCaugc]#N#g' hairpin.fa > hairpin_parse.fa
sed -i 's#\s.*##' hairpin_parse.fa
seqkit grep -r --pattern ".*hsa-.*" hairpin_parse.fa > hairpin_hsa.fa
seqkit seq --rna2dna hairpin_hsa.fa > tmp.fa
fasta_formatter -w 0 -i tmp.fa -o tmp1.fa
rm hairpin.fa hairpin_hsa.fa hairpin_parse.fa tmp.fa
mv tmp1.fa hairpin.fa

# Index miRBase hairpin file
bowtie-build hairpin.fa hairpin

# For each sample:
# 1. Trim adapters
# 2. Collapse reads
# 3. Align to hairpin (allow multi-hits to be reported up to 50 times)
# 4. Extract mapped reads to BAM
# 5. Extract unmapped reads to FASTQ
# 6. Repeat 2,3,4,5 until no more unmapped reads remain for sample.

while read -r line; do
    
    bn=$( basename $line .fastq.gz)
    sample=$( echo $bn | cut -d '-' -f1-2)
    
    mkdir ${sample}
    fastqc $line -o ${sample}
 
    trim_galore --adapter TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --length 17 --max_length 30 --three_prime_clip_R1 4 --clip_R1 4 --fastqc $line --output_dir ${sample}
    
    seqcluster collapse -f ${sample}/*trimmed.fq.gz -m 1 --min_size 15 -o collapsed 
    (bowtie -x hairpin -q collapsed/${sample}*trimmed.fastq -k 50 -e 99999 --best --strata --sam > ${sample}/${sample}_run1.sam) 2> ${sample}_run1.err

    processed=$(grep -o 'reads processed: [0-9]*' ${sample}_run1.err | awk '{print $3}') 
    alignment=$(grep -o 'reads with at least one alignment: [0-9]*' ${sample}_run1.err | awk '{print $7}') 
    failed=$(grep -o 'reads that failed to align: [0-9]*' ${sample}_run1.err | awk '{print $6}') 
    reported=$(grep -o 'Reported [0-9]* alignments' ${sample}_run1.err | awk '{print $2}')

    rm ${sample}_run1.err
    printf "%s %s %s %s %s %s\n" "$sample" "1" "$processed" "$alignment" "$failed" "$reported" > ${sample}_run1.err

    for ((i=1; i<=20; i++)); do
        
        samtools view -F 4 -h -b ${sample}/${sample}_run${i}.sam > ${sample}/${sample}_run${i}.bam 
        samtools sort -o ${sample}/${sample}_run${i}_srt.bam ${sample}/${sample}_run${i}.bam
        samtools view -f 4 ${sample}/${sample}_run${i}.sam | samtools fastq - > ${sample}/${sample}_run${i}_unmapped.fq
        output=$(samtools view -f 4 ${sample}/${sample}_run${i}.sam | samtools fastq - 2>&1)
    
            if [[ $output == *"processed 0 reads"* ]]; then
                echo "Finished"
                break
            else
                j=$((i + 1))
                seqcluster collapse -f ${sample}/${sample}_run${i}_unmapped.fq -m 1 --min_size 15 -o collapsed
                (bowtie -x hairpin -q collapsed/${sample}_run${i}*trimmed.fastq -k 50 -e 99999 --best --strata --sam --trim5 1 --trim3 1 > ${sample}/${sample}_run${j}.sam) 2> ${sample}_run${j}.err

                processed=$(grep -o 'reads processed: [0-9]*' ${sample}_run${j}.err | awk '{print $3}')
                alignment=$(grep -o 'reads with at least one alignment: [0-9]*' ${sample}_run${j}.err | awk '{print $7}') 
                failed=$(grep -o 'reads that failed to align: [0-9]*' ${sample}_run${j}.err | awk '{print $6}') 
                reported=$(grep -o 'Reported [0-9]* alignments' ${sample}_run${j}.err | awk '{print $2}')

                rm ${sample}_run${j}.err
                printf "%s %s %s %s %s %s\n" "$sample" "$j" "$processed" "$alignment" "$failed" "$reported" > ${sample}_run${j}.err
            fi
    done

    samtools merge ${sample}.merged.bam ${sample}/*_srt.bam
    samtools index ${sample}.merged.bam
    cat *run*.err > ${sample}.err && rm *run*.err
done < samples.txt

# Quantify the *merged.bam files for each sample:
mirtop gff --hairpin hairpin.fa --gtf hsa.gff -o mirtop --sps hsa *merged.bam
mirtop counts --hairpin hairpin.fa --gtf hsa.gff -o mirtop --sps hsa --add-extra --gff mirtop/mirtop.gff
mirtop export --format isomir --hairpin hairpin.fa --gtf hsa.gff --sps hsa -o mirtop mirtop/mirtop.gff
mirtop stats mirtop/mirtop.gff --out mirtop/stats
Rscript collapse_mirtop.r mirtop/mirtop.tsv
