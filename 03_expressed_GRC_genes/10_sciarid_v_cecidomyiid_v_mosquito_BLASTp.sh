#!/bin/bash -l

#SBATCH --job-name=cec_blastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --export=ALL
#SBATCH --time=0-02:00:00
#SBATCH --partition=ac3-compute
#SBATCH --mem=4gb
#SBATCH --output=cec_blastp.%j.log

set -e
hostname

source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics
 

WORKDIR=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/cecidomyiid_v_sciarid
mkdir -p ${WORKDIR}
cd ${WORKDIR}

# ── 1) Download outgroup ───────────────────────────────────────────────────────────────
# echo "Downloading Anopheles gambiae (AgamP4) RefSeq protein set..."
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/943/734/735/GCF_943734735.2_idAnoGambNW_F1_1/GCF_943734735.2_idAnoGambNW_F1_1_protein.faa.gz

# gunzip GCF_943734735.2_idAnoGambNW_F1_1_protein.faa.gz

OUTGROUP_AA=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/cecidomyiid_v_sciarid/GCF_943734735.2_idAnoGambNW_F1_1_protein.faa

# ── 2) Extract core-transcripts only from B_cop ─────────────────────────────────────────────────────────────────
BCOP_GFF=/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked_core_and_grc_braker3.gff3
BCOP_GENOME=/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked.fa

gffread \
    -g $BCOP_GENOME \
    -y Bcop_all_proteins.faa \
    $BCOP_GFF

seqkit grep -r -p "^s" Bcop_all_proteins.faa > Bcop_core_proteins.faa

SCI_AA=Bcop_core_proteins.faa
CEC_AA=/mnt/loki/ross/flies/cecidomyiidae/Aphidoletes_aphidimyza/aphid_midge_inversion/15_braker/outputs/M_producer_ragtag.Aphidoletes_aphidimyza.braker/Aaph_braker.aa

# ── 3) makde db ───────────────────────────────────────────────────────────────
echo "Building Cecidomyiid (Aphidoletes aphidimyza) DB..."
makeblastdb -in $CEC_AA -dbtype prot -parse_seqids -out Cecidomyiid_DB
 
echo "Building Sciarid core (GRC-excluded) DB..."
makeblastdb -in $SCI_AA -dbtype prot -parse_seqids -out Sciarid_core_DB
 
echo "Building Outgroup (Anopheles gambiae) DB..."
makeblastdb -in $OUTGROUP_AA -dbtype prot -parse_seqids -out Outgroup_DB
 
# ── 3) Three way BLASTp  ───────────────────────────────────────────────────────────────
SCRATCH=/scratch/${USER}/BLASTp_3way_GRC.${SLURM_JOB_ID}
mkdir -p ${SCRATCH}
cd ${SCRATCH}

echo "Copying input files..."
rsync -av ${WORKDIR}/GRC_proteins.faa ${SCRATCH}/      # 10 translated GRC queries
rsync -av ${WORKDIR}/Cecidomyiid_DB.*   ${SCRATCH}/
rsync -av ${WORKDIR}/Sciarid_core_DB.*  ${SCRATCH}/
rsync -av ${WORKDIR}/Outgroup_DB.*      ${SCRATCH}/
 
run_blast () {
  local db=$1
  local out=$2
  echo "Running BLASTp vs ${db}..."
  blastp \
    -query GRC_proteins.faa \
    -db ${db} \
    -out ${out} \
    -evalue 10 \
    -num_threads 16 \
    -max_target_seqs 10 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'
}
 
run_blast Cecidomyiid_DB  GRC_vs_Cecidomyiid.tsv
run_blast Sciarid_core_DB GRC_vs_SciaridCore.tsv
run_blast Outgroup_DB     GRC_vs_Outgroup.tsv
 
echo "Syncing results back..."
rsync -av *.tsv ${WORKDIR}/
 
echo "Cleaning scratch directory..."
rm -rf ${SCRATCH}

# ── 4) get Similarity Index  ───────────────────────────────────────────────────────────────
cd ${WORKDIR}
SCRIPT=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/scripts/alien_index.py

python $SCRIPT \
        --cecidomyiid GRC_vs_Cecidomyiid.tsv \
        --sciarid     GRC_vs_SciaridCore.tsv \
        --outgroup    GRC_vs_Outgroup.tsv \
        --out         Alien_Index_summary.csv

echo "Done."
