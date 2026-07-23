import os
import argparse
from Bio import SeqIO

"""
Calculate GC content % in W bp window using sliding step n
"""

def gc_content(seq):
    gc = seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0

def at_content(seq):
    at = seq.count("A") + seq.count("T") + seq.count("a") + seq.count("t")
    return (at / len(seq)) * 100 if len(seq) > 0 else 0

def compute_gc_sliding(seq, window_size, step):
    pos_gc_dict = {}
    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        gc = gc_content(window)
        at = at_content(window)
        position = int(i + window_size // 2)  # midpoint of window
        pos_gc_dict[position] = {"GC":gc, "AT":at}
    return pos_gc_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate GC content in n bp window using sliding window of size w")
    parser.add_argument("-g", "--genome", type=str, help="Path to genome fasta file", required=True)
    parser.add_argument("-w", "--window", type=str, help="Window_size", required=False, default=2000)
    parser.add_argument("-n", "--step_size", type=str, help="Step size", required=False, default=1000)
    args = parser.parse_args()

    window_size = int(args.window)
    step = int(args.step_size)

    #get filename for output
    basename = os.path.basename(args.genome)
    # Parse genome
    output_file = f"GC_{basename}_w{window_size}_n{step}.tsv"
    with open(output_file, "w") as outf:
        outf.write("Scaffold\tPosition\tGC_Content\tAT_Content\n")
        for record in SeqIO.parse(args.genome, "fasta"):
            print(f"Counting GC content for {record.id} in window size {window_size} with step {step} for {len(record.seq)} bp")
            pos_gc_dict = compute_gc_sliding(str(record.seq), window_size, step)
            for pos, percentage in pos_gc_dict.items():
                gc = percentage["GC"]
                at = percentage["AT"]
                outf.write(f"{record.id}\t{pos}\t{gc:.2f}\t{at:.2f}\n")
