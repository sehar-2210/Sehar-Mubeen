from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO

# Alignment function
def align_seq(seq1, seq2, match=2, mismatch=-2, gap_open=-2, gap_extend=-1):
    alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_extend)
    best_alignment = alignments[0]
    
    print("The alignment score is: ", best_alignment.score)
    print("The alignment sequence 1 is: ", best_alignment.seqA)
    print("The alignment sequence 2 is: ", best_alignment.seqB)
    print("The start of alignment is: ", best_alignment.start)
    print("The end of alignment is: ", best_alignment.end)

    return best_alignment

# Similarity function
def similarity(alignment):
    seq1 = alignment.seqA
    seq2 = alignment.seqB
    start = alignment.start
    end = alignment.end
    aligned1 = seq1[start:end]
    aligned2 = seq2[start:end]

    matches = 0
    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i] and aligned1[i] != '-':
            matches += 1

    print("matches: ", matches)
    length = end - start
    similarity_score = (matches / length) * 100 if length > 0 else 0
    print("The similarity of the alignment is:", round(similarity_score, 2), '%')
    return similarity_score

# Gap frequency function
def gap_frequency(alignment):
    seq1 = alignment.seqA
    seq2 = alignment.seqB

    gaps_seq1 = seq1.count('-')
    gaps_seq2 = seq2.count('-')
    total_gaps = gaps_seq1 + gaps_seq2
    alignment_length = len(seq1)

    gap_freq = (total_gaps / alignment_length) * 100 if alignment_length > 0 else 0

    print("Gaps in sequence 1:", gaps_seq1)
    print("Gaps in sequence 2:", gaps_seq2)
    print("Total gaps:", total_gaps)
    print("Gap frequency:", round(gap_freq, 2), '%')
    return gap_freq

# Conserved regions (≥ 20bp)
def find_conserved_regions(alignment, threshold=20):
    seq1 = alignment.seqA
    seq2 = alignment.seqB

    conserved = []
    current_region = ""
    start_index = None

    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a == b and a != "-":
            if start_index is None:
                start_index = i
            current_region += a
        else:
            if len(current_region) >= threshold:
                conserved.append((start_index, start_index + len(current_region), current_region))
            current_region = ""
            start_index = None

    if len(current_region) >= threshold:
        conserved.append((start_index, start_index + len(current_region), current_region))

    print(f"\n🔍 Conserved Regions (length ≥ {threshold} bp):")
    if conserved:
        for i, (start, end, region) in enumerate(conserved, 1):
            print(f"{i}. {region} (start: {start}, end: {end}, length: {len(region)})")
    else:
        print("No conserved regions found.")

    return conserved

# ✅ Read sequences from FASTA files
def read_fasta_sequence(file_path):
    record = SeqIO.read(file_path, "fasta")
    return record.seq

# ✅ MAIN
if __name__ == "__main__":
    seq1 = read_fasta_sequence("seq1.fasta")
    seq2 = read_fasta_sequence("seq2.fasta")

    alignment_results = align_seq(seq1, seq2)
    similarity(alignment_results)
    gap_frequency(alignment_results)
    find_conserved_regions(alignment_results, threshold=20)
