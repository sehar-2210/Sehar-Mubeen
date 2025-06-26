import csv

def is_valid_dna(sequence):
    valid_nucleotides = set("ATGC")
    return set(sequence.upper()).issubset(valid_nucleotides)

def calculate_gc_content(sequence):
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return round(gc_count / len(sequence) * 100, 2) if sequence else 0.0

def read_fasta(filename):
    sequences = {}
    with open(filename, 'r') as file:
        current_id = ""
        current_seq = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    return sequences

def find_unique_nucleotides(sequences_dict):
    all_nucleotides = set()
    for seq in sequences_dict.values():
        all_nucleotides.update(seq.upper())
    return all_nucleotides

def analyze_and_save(fasta_file, output_csv):
    sequences = read_fasta(fasta_file)

    if not sequences:
        print("⚠ No sequences found in the FASTA file.")
        return

    print(f"✅ Total sequences found: {len(sequences)}")
    unique_nucleotides = find_unique_nucleotides(sequences)
    print(f"🔬 Unique nucleotides found: {unique_nucleotides}")

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ID", "Length", "GC_Content(%)", "Is_Valid"])

        for seq_id, seq in sequences.items():
            gc_content = calculate_gc_content(seq)
            valid = is_valid_dna(seq)
            writer.writerow([seq_id, len(seq), gc_content, valid])
            print(f"📄 Processed: {seq_id} | Length: {len(seq)} | GC: {gc_content}% | Valid: {valid}")

    print(f"\n✅ Analysis complete! CSV saved at: {output_csv}")

if __name__ == "__main__":
    fasta_input = "sequence.fasta"
    csv_output = "sequence_analysis.csv"
    analyze_and_save(fasta_input, csv_output)
