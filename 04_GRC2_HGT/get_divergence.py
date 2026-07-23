'''
This script is used to estimate the age of the HGT region in GRC
 
We use FastGA (https://github.com/thegenemyers/FASTGA?tab=readme-ov-file#fastga-a-fast-genome-aligner) 
to align the Rickettsiaceae genome to GRC2 and claculate the differences between the alignmnets (df:i)

Given divergence time (T) between 2 species:

T = D / (μ1 + μ2)

where 
    μ1 = mutation rate in species 1
    μ2 = mutation rate in species 2


We can sum the divergence for each alignment block and get a crude estimate for T 

'''
# Usage : python3 get_divergence.py

# --- mutation rates ---
MU_FLY = 2.8e-9  # per site per generation 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC3872194/
MU_BACTERIA = 1.1e-9  # per site per generation
# https://academic.oup.com/gbe/article/10/3/723/4838064?login=false

PAF_FILE = "grc2_vs_rickettsia_1to1.1aln.paf"


def get_divergence(alignment_file):
    """Compute divergence from df:i (counted differences)"""
    total_df = 0
    total_len = 0

    for line in alignment_file:
        columns = line.strip().split('\t')
        
        # get alignment length from column 10
        alignment_length = int(columns[9])
        # parse df:i
        differences = int(columns[13].split(":")[2])

        total_df += differences
        total_len += alignment_length

    return total_df / total_len

def estimate_hgt_age(divergence, mu_fly=MU_FLY, mu_bacteria=MU_BACTERIA):
    """Estimate HGT age in years"""
    return divergence / (mu_fly + mu_bacteria)


def main():
    with open(PAF_FILE, "r") as f:
        D = get_divergence(f)
    
    T = estimate_hgt_age(D)

    print(f"Divergence: {D:.6f}")
    print(f"HGT age estimate: {T:.0f} generations")

if __name__ == "__main__":
    main()
