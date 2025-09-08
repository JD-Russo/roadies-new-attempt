from collections import defaultdict

def remove_duplicate_fasta_headers(input_file, output_file):
    seen_headers = set()
    duplicates = defaultdict(int)
    kept_sequences = []

    with open(input_file, 'r') as infile:
        current_header = None
        current_seq = []

        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header:
                    if current_header not in seen_headers:
                        seen_headers.add(current_header)
                        kept_sequences.append((current_header, ''.join(current_seq)))
                    else:
                        duplicates[current_header] += 1

                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Handle the last sequence
        if current_header:
            if current_header not in seen_headers:
                kept_sequences.append((current_header, ''.join(current_seq)))
            else:
                duplicates[current_header] += 1

    # Write unique sequences to output
    with open(output_file, 'w') as outfile:
        for header, sequence in kept_sequences:
            outfile.write(f"{header}\n")
            # Split sequence into lines of 60 characters
            for i in range(0, len(sequence), 60):
                outfile.write(sequence[i:i+60] + "\n")

    print(f"Found {sum(duplicates.values())} duplicate headers.")
    print(f"{len(duplicates)} unique headers had duplicates.")
    print(f"Cleaned file written to: {output_file}")

# Example usage
remove_duplicate_fasta_headers("out.fa", "cleaned_out.fa")
