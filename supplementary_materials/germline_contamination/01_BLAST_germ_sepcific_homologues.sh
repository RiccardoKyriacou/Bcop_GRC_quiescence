#!/bin/bash -l

#SBATCH --job-name=BLAST
#SBATCH --nodes=1 #number of nodes requested
#SBATCH --ntasks=16 #number of threads per node
#SBATCH --export=ALL
#              d-hh:mm:ss
#SBATCH --time=0-90:00:00 # Upper time limit for job
#SBATCH --partition ac3-compute
#SBATCH --mem=32gb # How much memory you need.
#SBATCH --output=BLAST.%j.log   # Standard output and error log
#SBATCH --mail-type=END,FAIL # Turn on mail notification
#SBATCH --mail-user=s2673271@ed.ac.uk  # Where to send mail

# Exit script on error
set -e

# Define and create a unique scratch directory for this job
SCRATCH=/scratch/${USER}/BLAST.${SLURM_JOB_ID}
mkdir -p $SCRATCH
cd $SCRATCH

# Activate the conda environment for the job
source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

# Sync in files 
# Gemrline genes to query 
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/05_contamiantion/outputs/germline_CDS.fasta $SCRATCH
# Core and GRC fasta files to make blast db 
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/05_contamiantion/outputs/annotated_genes.fasta $SCRATCH


### STEP 1) BLAST GRC genes against core genes ###

# Make a blast data base of all annotated genes 
echo "making BLASTDB"
makeblastdb \
-in annotated_genes.fasta \
-dbtype nucl \
-parse_seqids \
-out annotated_genes_DB

echo ">ls"
ls 

# Blast GRC genes against core genes 
echo "BLASTn..."
tblastx \
-query germline_CDS.fasta \
-db annotated_genes_DB \
-out germ_specific_homologue_blast_output.tsv \
-outfmt '6 std qlen slen qseq sseq'

# Sync outputs out
echo "Syncing results back..."
rsync -av $SCRATCH/*_output.tsv /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/05_contamiantion/outputs/

# Clear and delete scratch
rm -rf ${SCRATCH}

# Finish the script
exit 0
