"""Call natural-library PETase blocks from ONT FASTQ.gz reads.

Adapted from PETase_ONT.zip. Mapping thresholds, segment boundaries, parent
tie-breaking, and output columns follow the original pipeline.
"""
import argparse
import gzip
from pathlib import Path
import subprocess
import tempfile

from Bio import Align, SeqIO
from Bio.Seq import Seq

# Biological blocks alternate with Golden Gate junction segments.
BOUNDARIES = [0, 159, 183, 244, 271, 333, 357, 408, 432, 504,
              528, 589, 613, 676, 700, 769, 793, 850, 874, 1046]


def read_alignment(path):
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    return [''.join(record.splitlines()[1:]).upper()
            for record in '\n'.join(lines).split('>') if record]


def pairwise_aligner(extension):
    # Explicitly preserve historical Biopython defaults: match=1, mismatch=0.
    # The original local variable mismatch_score=-1 never changed the aligner.
    aligner = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=0,
                                  open_gap_score=-10, extend_gap_score=extension)
    aligner.end_gap_score = 0
    return aligner


def map_reads(fastq, reference, output):
    parents = read_alignment(reference)
    if len(parents) != 10 or {len(p) for p in parents} != {1047}:
        raise ValueError('Expected ten aligned parents of length 1047.')
    consensus = ''.join(max('ACGT-', key=column.count) for column in zip(*parents))
    orienter = pairwise_aligner(-0.1)
    quality_aligner = pairwise_aligner(-1)
    processed = called = 0

    # Each invocation gets its own temporary files, so jobs cannot overwrite
    # another job's alignment inputs or outputs.
    with tempfile.TemporaryDirectory(prefix='petase_ont_') as directory:
        scratch = Path(directory)
        parent_file = scratch / 'parents.fasta'
        parent_file.write_text(''.join(f'>seq{i}\n{p}\n' for i, p in enumerate(parents)))
        read_file = scratch / 'read.fasta'
        aligned_file = scratch / 'aligned.fasta'
        with gzip.open(fastq, 'rt') as handle, Path(output).open('w') as outfile:
            for record in SeqIO.parse(handle, 'fastq'):
                processed += 1
                read = str(record.seq)
                if len(read) <= 500:
                    continue
                reverse = str(Seq(read).reverse_complement())
                forward_score = orienter.score(consensus, read)
                reverse_score = orienter.score(consensus, reverse)
                if forward_score > reverse_score and forward_score > 350:
                    pass
                elif reverse_score > forward_score and reverse_score > 350:
                    read = reverse
                else:
                    continue

                read_file.write_text(f'>NGSread\n{read}\n')
                with aligned_file.open('w') as out:
                    subprocess.run(['mafft', '--quiet', '--thread', '1', '--add',
                                    str(read_file), '--keeplength', str(parent_file)],
                                   stdout=out, check=True)
                alignment = read_alignment(aligned_file)
                if len(alignment) != 11 or len(alignment[-1]) != len(parents[0]):
                    raise ValueError('Unexpected MAFFT alignment shape.')
                calls = []
                for left, right in zip(BOUNDARIES, BOUNDARIES[1:]):
                    fragment = alignment[-1][left:right]
                    scores = [sum(a == b for a, b in zip(p[left:right], fragment))
                              for p in alignment[:-1]] + [fragment.count('-')]
                    calls.append(scores.index(max(scores)))

                chimera = ''.join(parents[parent][left:right] if parent < 10
                                  else '-' * (right-left)
                                  for parent, left, right in
                                  zip(calls, BOUNDARIES, BOUNDARIES[1:])).replace('-', '')
                matched = quality_aligner.align(chimera, read)[0]
                mismatch = matched.counts().mismatches
                gaps = ''.join('-' if p == -1 else 'X' for p in matched.indices.min(0))
                indels = gaps.strip('-').count('-')
                call_string = ''.join(str(parent) if parent < 10 else '-' for parent in calls)
                outfile.write(f'@{record.id},{call_string},{mismatch},{indels}\n')
                called += 1
    print(f'{processed} input reads; {called} mapped reads; output: {output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fastq', type=Path, help='One demultiplexed natural-library FASTQ.gz file')
    parser.add_argument('--reference', type=Path,
                        default=Path(__file__).with_name('natural_aligned.fasta'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not args.fastq.name.endswith('.fastq.gz'):
        parser.error('Expected a .fastq.gz input file.')
    output = args.output or Path(args.fastq.name[:-9] + '_parsed_data.txt')
    map_reads(args.fastq, args.reference, output)
