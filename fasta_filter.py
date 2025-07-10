# Function to read FASTA file and validate content
def read_fasta(file_path):
    sequences = {}
    try:
        with open(file_path, 'r') as file:
            header = ''
            seq = ''
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header and seq:
                        sequences[header] = seq
                    header = line
                    seq = ''
                else:
                    if not all(c in 'ATGCatgc' for c in line):
                        raise ValueError(f" Invalid characters in line: {line}")
                    seq += line
            if header and seq:
                sequences[header] = seq
        return sequences
    except FileNotFoundError:
        print(f" File not found: {file_path}")
        return {}
    except Exception as e:
        print(f" Error reading FASTA: {e}")
        return {}

# Function to filter sequences by minimum length
def filter_sequences(sequences, min_length):
    filtered = {h: s for h, s in sequences.items() if len(s) >= min_length}
    return filtered

# Function to write filtered sequence
def write_fasta(sequences, output_path):
    try:
        with open(output_path, 'w') as out:
            for header, seq in sequences.items():
                out.write(f"{header}\n")
                for i in range(0, len(seq), 60):
                    out.write(f"{seq[i:i+60]}\n")
    except Exception as e:
        print(f" Error writing to file: {e}")
def main():
    input_file = input(" Enter input FASTA file name (e.g., sample.fasta): ").strip()
    sequences = read_fasta(input_file)

    if not sequences:
        return

    total_read = len(sequences)

    try:
        min_len = int(input("📏 Enter minimum sequence length to keep: "))
    except ValueError:
        print(" Invalid number entered.")
        return

    filtered = filter_sequences(sequences, min_len)
    total_passed = len(filtered)

    output_file = input(" Enter output FASTA file name (e.g., filtered.fasta): ").strip()
    write_fasta(filtered, output_file)

    print("\n Summary Statistics:")
    print(f" Total sequences read: {total_read}")
    print(f" Sequences passing filter: {total_passed}")

if __name__ == "__main__":
    main()
