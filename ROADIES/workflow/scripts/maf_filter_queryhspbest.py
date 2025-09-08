#!/usr/bin/env python3
"""maf_filter_queryhspbest.py

Filters a MAF file so that for each query gene, only the top-N HSPs (by score)
are kept. Any HSPs scoring below the Nth best for that query are discarded.

Usage:
  python maf_filter_queryhspbest.py input.maf --queryhspbest 20 --output filtered.maf
"""
import argparse
import gzip
import pathlib
import re
import sys
from collections import defaultdict
from typing import List, Tuple


def open_maybe_gzip(path: pathlib.Path):
    return gzip.open(path, 'rt') if path.suffix == '.gz' else open(path, 'rt')


def parse_maf_blocks(path: pathlib.Path) -> List[Tuple[int, str, str]]:
    """
    Parse the input MAF into a list of tuples:
      (block_score, query_gene_id, raw_block_text)
    where raw_block_text is the full text of the block, including the 'a ' line
    and all 's ' lines, ending with a blank line.
    """
    blocks = []
    gene_id_pat = re.compile(r"gene_[0-9]+", re.I)
    score_pat = re.compile(r"score=([-0-9]+)")

    with open_maybe_gzip(path) as fh:
        buffer_lines = []
        current_score = None
        query_gene = None
        collecting = False

        for line in fh:
            buffer_lines.append(line)
            stripped = line.rstrip("\n")

            # Detect 'a score=' line
            if stripped.startswith('a '):
                m = score_pat.search(stripped)
                current_score = int(m.group(1)) if m else 0
                collecting = True
                continue

            # If we are inside a block, look for the first 's ' that gives query gene
            if collecting and stripped.startswith('s '):
                parts = stripped.split(maxsplit=6)
                src = parts[1]
                m_gid = gene_id_pat.search(src)
                if m_gid and query_gene is None:
                    query_gene = m_gid.group(0)
                continue

            # Blank line indicates end of block
            if stripped == '':
                if collecting and current_score is not None and query_gene:
                    raw_text = ''.join(buffer_lines)
                    blocks.append((current_score, query_gene, raw_text))
                # Reset for next block
                buffer_lines = []
                current_score = None
                query_gene = None
                collecting = False

        # Handle last block if file didn't end with blank
        if collecting and current_score is not None and query_gene:
            raw_text = ''.join(buffer_lines)
            blocks.append((current_score, query_gene, raw_text))

    return blocks


def filter_top_n(blocks: List[Tuple[int, str, str]], top_n: int) -> List[str]:
    """
    Given a list of (score, query_gene, raw_block_text), filter to keep only
    those blocks whose score is among the top-N for that query_gene.

    Return a list of raw_block_text strings to keep.
    """
    # Group scores by query_gene
    scores_by_gene = defaultdict(list)
    for score, gene, _ in blocks:
        scores_by_gene[gene].append(score)

    # Determine threshold per gene: the Nth largest score (if fewer than N, keep all)
    threshold_by_gene = {}
    for gene, scores in scores_by_gene.items():
        sorted_scores = sorted(scores, reverse=True)
        if len(sorted_scores) < top_n:
            threshold = sorted_scores[-1]
        else:
            threshold = sorted_scores[top_n - 1]
        threshold_by_gene[gene] = threshold

    # Now select blocks whose score >= threshold for their gene
    kept_blocks = []
    for score, gene, raw_text in blocks:
        if score >= threshold_by_gene[gene]:
            kept_blocks.append(raw_text)
    return kept_blocks


def write_filtered_maf(block_texts: List[str], out_path: pathlib.Path):
    with open(out_path, 'w') as fh:
        for block in block_texts:
            fh.write(block)
            if not block.endswith("\n"):
                fh.write("\n")


def build_cli():
    parser = argparse.ArgumentParser(description="Filter MAF by queryHSPbest (top-N per query gene).")
    parser.add_argument('input', type=pathlib.Path,
                        help='Input MAF file (possibly .gz)')
    parser.add_argument('--queryhspbest', type=int, default=20,
                        help='Retain only top-N HSPs per query gene')
    parser.add_argument('--output', type=pathlib.Path, required=True,
                        help='Path to write filtered MAF')
    return parser


def main(argv=None):
    parser = build_cli()
    args = parser.parse_args(argv)

    blocks = parse_maf_blocks(args.input)
    kept = filter_top_n(blocks, args.queryhspbest)
    write_filtered_maf(kept, args.output)
    print(f"Filtered MAF: kept {len(kept)} blocks (out of {len(blocks)})")

if __name__ == '__main__':
    sys.exit(main())

