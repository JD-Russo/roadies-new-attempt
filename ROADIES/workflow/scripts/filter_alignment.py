import argparse
import multiprocessing as mp

# Fixed thresholds
IDENTITY_THRESHOLD = 65.0
COVERAGE_THRESHOLD = 85.0
CONTINUITY_THRESHOLD = 85.0

def parse_maf_block(block):
    """Extract sequences and sizes from a maf block."""
    lines = [line for line in block.strip().split('\n') if line.startswith('s')]
    if len(lines) != 2:
        return None  # Only handle pairwise alignments

    seqs = []
    src_sizes = []
    aligned_sizes = []

    for line in lines:
        parts = line.split()
        src = parts[1]
        start = int(parts[2])
        size = int(parts[3])
        strand = parts[4]
        srcSize = int(parts[5])
        text = parts[6]
        seqs.append(text)
        src_sizes.append(srcSize)
        aligned_sizes.append(size)

    return seqs, src_sizes, aligned_sizes

def compute_identity(seq1, seq2):
    """Compute identity = matches / non-gap columns."""
    matches = 0
    aligned_positions = 0
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            aligned_positions += 1
            if a.casefold() == b.casefold():
                matches += 1
    if aligned_positions == 0:
        return 0.0
    return matches / aligned_positions * 100

def compute_continuity(seq1, seq2):
    """Compute continuity = non-gap columns / total columns."""
    non_gap_columns = 0
    total_columns = len(seq1)
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            non_gap_columns += 1
    if total_columns == 0:
        return 0.0
    return non_gap_columns / total_columns * 100

def compute_coverage(aligned_sizes, src_sizes):
    """Compute coverage = aligned size / shorter total size."""
    aligned_length = min(aligned_sizes)
    shorter_seq_len = min(src_sizes)
    if shorter_seq_len == 0:
        return 0.0
    return aligned_length / shorter_seq_len * 100

def process_block(block):
    """Process a single maf block and return it if it passes filters."""
    if not block.strip():
        return None

    parsed = parse_maf_block(block)
    if parsed is None:
        return None

    seqs, src_sizes, aligned_sizes = parsed

    identity = compute_identity(seqs[0], seqs[1])
    continuity = compute_continuity(seqs[0], seqs[1])
    coverage = compute_coverage(aligned_sizes, src_sizes)

    if (identity >= IDENTITY_THRESHOLD and
        continuity >= CONTINUITY_THRESHOLD and
        coverage >= COVERAGE_THRESHOLD):
        return block.strip()
    else:
        return None

def filter_maf_parallel(input_file, output_file, num_workers=None):
    """Filter maf file using multiprocessing."""
    with open(input_file) as f:
        content = f.read()

    blocks = content.strip().split('\n\n')
    print(f"Total alignment blocks: {len(blocks)}")

    # Create a multiprocessing Pool
    with mp.Pool(processes=num_workers) as pool:
        results = pool.map(process_block, blocks)

    # Write out only non-None results
    with open(output_file, 'w') as out:
        for block in results:
            if block:
                out.write(block + '\n\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter MAF file with multiprocessing.")
    parser.add_argument("--input", required=True, help="Input MAF file")
    parser.add_argument("--output", required=True, help="Output filtered MAF file")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes (default: 4)")

    args = parser.parse_args()

    filter_maf_parallel(args.input, args.output, num_workers=args.threads)