# Megascale recombination of PETase

Construct preparation and publication analyses for the natural-parent PETase recombination library. The three notebooks are self-contained and include their executed outputs.

## Notebooks

| Notebook | What it does | Outputs |
| --- | --- | --- |
| [design_DNA_constructs.ipynb](design_DNA_constructs.ipynb) | Recreates the ten Twist synthesis constructs, including Golden Gate junctions and vector sequences | `constructs/*.fasta` |
| [analyze_sequencing.ipynb](analyze_sequencing.ipynb) | Explores sequencing cutoffs, computes Figure 3 enrichment and diversity, and plots the saved Figure 4C identity matrix | `figures/N488e.eps`, `NPEe.eps`, `NPE_pwdist_hist.eps`, `PW_iden.eps` |
| [analyze_activity.ipynb](analyze_activity.ipynb) | Recreates the screening and reaction-progress panels from the original pasted tables | `figures/screening_data.eps`, `reaction_progress.eps` |

These follow the original `Twist_designs_final.ipynb`, `paper_analysis.ipynb`, and `process_pete_data.ipynb`, respectively. Comments explain the retained analysis choices.

## Run

Use Python 3.12 and start in the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Open a notebook in your Jupyter frontend or notebook editor, select this environment, and run all cells. The notebooks are independent. They regenerate the FASTA or EPS outputs and display the plots inline. Reading the rendered notebooks on GitHub requires no installation.

## Files

- `data/constructs/`: the ten fixed-breakpoint parent rows, their saved Twist codon-optimization results, and nine Golden Gate junctions. Parent row 0 corresponds to synthesis record `natural1`.
- `data/sequencing/`: unchanged final SI block table and the three original sample-call files, losslessly gzip-compressed. The analysis reads these directly.
- `data/assays/`: original Excel workbooks, with a short note identifying the source ranges and rounding used in the plots.
- `assets/ssn/`: original RepNode50 XGMML network, candidate image `PETase80_5.png`, and transcribed run parameters. Open the XGMML in Cytoscape.
- `assets/structure/petase_recombination_blocks.pse`: the original PyMOL session, renamed from `structure.pse`. Open it in PyMOL.
- `ngs/`: natural-library ONT mapping script, aligned reference, and local/HTCondor running instructions.


## License

[MIT License](LICENSE). Copyright (c) 2026 RomeroLab.
