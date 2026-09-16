# ONT read mapping

`fastq2chimera.py` is adapted from the recovered `PETase_ONT.zip` pipeline for the natural-parent library. It reads a demultiplexed FASTQ.gz file and maps each accepted read to the parent blocks using the original aligned reference and MAFFT.

## Local run

Install the repository's Python requirements and MAFFT, available as `mafft` on PATH. Then, from the repository root:

```bash
mkdir -p ngs/output
python ngs/fastq2chimera.py /path/to/barcode11/read_file.fastq.gz \
  --output ngs/output/read_file_parsed_data.txt
```

The same command applies to each natural-library FASTQ from barcodes 11–13. The mapping reference defaults to `natural_aligned.fasta` beside the script; it is an **input**, not an output. It includes the assembly junctions and flanks, so the SI block table used for downstream sequence reconstruction is not a drop-in replacement.

The script keeps reads longer than 500 nt and selects the uniquely better orientation if its pairwise alignment score exceeds 350. MAFFT adds the oriented read to the ten aligned parents with `--add --keeplength`. At each segment, the highest-identity parent is chosen; an all-gap candidate can instead produce `-`. Ties use the first parent, as in the original code.

## Output and downstream analysis

One headerless CSV row is written per mapped read:

```text
read_id,19_segment_calls,mismatch_count,internal_indel_count
```

The 19 one-character calls alternate between ten biological blocks and nine junction segments. Biological blocks occupy zero-based positions 0, 2, …, 18. `-` denotes a missing segment and may occur at a biological block; it must not simply be deleted. Mismatches and internal gap bases are reported, not filtered by the mapper.

The archive contains mapping and job-submission code, but no verified filtering/aggregation script that produces `sample_NST.txt`, `sample_N488.txt`, and `sample_NPE.txt`. The exact quality thresholds and barcode-to-sample assignment therefore remain to be recorded. The analysis notebook starts from the supplied sample files and does not invent those steps.

## HTCondor

`run.sh` and `run.sub` retain the original one-file-per-job approach with simpler environment handling. Workers must already provide Python, Biopython, and MAFFT via PATH (or the cluster's standard environment). Submit from `ngs/`:

```bash
mkdir -p reads logs output
# Put only the natural-library FASTQ.gz files in reads/.
find reads -name '*.fastq.gz' -exec basename '{}' \; > inputs.txt
condor_submit run.sub
```

Input basenames must be unique. HTCondor transfers each read file plus the script and reference; resulting `*_parsed_data.txt` files return to the submission directory. Local mapping is tested; this submission template has not been run on the original cluster. The archived cluster used separately built Python/package bundles, which were not included in the ZIP and are not required for the local command.

## Changes and validation

The original mapping thresholds, 19 segment boundaries, parent tie-breaking, and output fields are retained. The reference has 1,047 aligned columns; the original boundary list ends at 1,046 and is preserved. File handling now uses arguments and per-job temporary directories, checks MAFFT's exit status, and streams FASTQ records. Pairwise scoring is explicit: match 1, mismatch 0, gap-open −10, gap-extension −0.1 for orientation or −1 for read-to-chimera comparison, with free terminal gaps. These preserve the effective original settings; its standalone `mismatch_score = -1` variable did not change the aligner.

Three natural reads reproduce the original saved calls, mismatch counts, and indel counts exactly using Python 3.12, Biopython 1.88, and MAFFT 7.526. This is a smoke test, not a full sequencing rerun. The historical environment instructions name Python 3.9 and MAFFT 7.505; an exact historical Biopython version was not recorded.
