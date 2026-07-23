#!/bin/bash -l

#SBATCH --job-name=STAR_TPM
#SBATCH --nodes=1 #number of nodes requested
#SBATCH --ntasks=16 #number of threads per node
#SBATCH --export=ALL
#              d-hh:mm:ss
#SBATCH --time=00-120:00:00 # Upper time limit for job
#SBATCH --partition ac3-compute
#SBATCH --mem=32gb # How much memory you need.
#SBATCH --output=STAR_TPM.%j.log   # Standard output and error log
#SBATCH --mail-type=END,FAIL # Turn on mail notification
#SBATCH --mail-user=s2673271@ed.ac.uk  # Where to send mail

# Exit script on error
set -e

# Define and create a unique scratch directory for this job
SCRATCH=/scratch/${USER}/STAR.${SLURM_JOB_ID}
mkdir -p $SCRATCH
cd $SCRATCH

# Activate the conda environment for the job
source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

# -------------------------------------------------------------------
# STEP 1: Sync files in 
# -------------------------------------------------------------------

# GRC genome
rsync -av /mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked.fa $SCRATCH
# Genome annotation (complete genome annotation)
rsync -av /mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked_core_and_grc_braker3.gff3 $SCRATCH
# Python script to parse StringTie output
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/01_STAR_Relaxed_Nmax2/scripts/get_TPM_values.py $SCRATCH

# copy embryo dir and unzip
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/Urban_RNA-seq_data/embryo_fastqs/*.fq.gz $SCRATCH

echo "---------------------------------------------------------------"
echo "Listing contents of $SCRATCH"
ls $SCRATCH
echo "---------------------------------------------------------------"

# -------------------------------------------------------------------
# STEP 2: Trim Reads
# -------------------------------------------------------------------
# fastp to trim  all reads
echo "Trimming reads with fastp..."
for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
	echo "Trimming $base"
  	fastp -i ${base}_1.fq.gz -I ${base}_2.fq.gz -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz && \
  	rm ${base}_1.fq.gz ${base}_2.fq.gz
done

# -------------------------------------------------------------------
# STEP 3: Index B coprophila genome
# -------------------------------------------------------------------

# Index the B_coprophila genome 
echo "Indexing B_coprophila genome..."
echo "run genomeGenerate"
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeSAindexNbases 12 \
--outFileNamePrefix idBraCopr2.1.primary.masked.fa \
--sjdbGTFfile idBraCopr2.1.primary.masked_core_and_grc_braker3.gff3  \
--genomeDir idBraCopr2.1.primary.masked.fa.STAR \
--genomeFastaFiles idBraCopr2.1.primary.masked.fa

# -------------------------------------------------------------------
# STEP 4: Map RNAseq reads to indexed B coprophila genome
#         --outFilterMismatchNmax 2
# -------------------------------------------------------------------

# Map RNAseq reads to B_coprophila genome
echo "Aligning B_coprophila RNAseq reads..."
for file in $(ls *_1.trimmed.fq.gz)
do
    base=$(basename "$file" "_1.trimmed.fq.gz")
    echo "Aligning $base to B_coprophila genome"
    STAR \
        --runThreadN 16 \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --readFilesIn ${base}_1.trimmed.fq.gz ${base}_2.trimmed.fq.gz \
        --outTmpDir ${base}.out \
        --outFileNamePrefix ${base}.STAR. \
        --outFilterMismatchNmax 2 \
        --genomeDir idBraCopr2.1.primary.masked.fa.STAR
done

# -------------------------------------------------------------------
# STEP 5: Post-mapping filtering for uniquely mapped reads
# -------------------------------------------------------------------

# Filter for any reads with a MAPQ=255 and are therefore are uniquely mapped
echo "Samtools MAPQ=255" 
for file in $(ls *.bam)
do
    base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
    output_file=${base}_uniquely_mapped.bam
    echo "Filtering unique reads for $file"
    samtools view -@ 16 -q 255 -b $file > $output_file
done

samtools index -M *_uniquely_mapped.bam

# -------------------------------------------------------------------
# STEP 6: Generate TPM estimates with StringTie
# -------------------------------------------------------------------

# stringtie TPM estimates
echo "stringtie TPM estimates"
for file in $(ls *_uniquely_mapped.bam)
do
    base=$(basename $file "_uniquely_mapped.bam")
    output_file=${base}_unique_mapping_TPM.gtf
    echo "Calculating TPM for $file"
    stringtie -o $output_file -G idBraCopr2.1.primary.masked_core_and_grc_braker3.gff3  $file
done

# -------------------------------------------------------------------
# STEP 7: Parse StringTie output using custom py script
# -------------------------------------------------------------------
python3 get_TPM_values.py -t .

# -------------------------------------------------------------------
# STEP 8: Re-sync output and delete scratch 
# -------------------------------------------------------------------

echo "Syncing mapped files and results back..."
rsync -av \
    *_uniquely_mapped.bam* \
    *.gtf \
    *.STAR.Log.final.out* \
    *.tsv \
    /mnt/loki/ross/flies/sciaridae/GRCs/Urban_RNA-seq_data


# Remove input files
echo "Removing input files"

# Clear and delete scratch
rm -rf ${SCRATCH}