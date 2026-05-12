# IUPred2A Nextflow

Nextflow wrapper for the bundled IUPred2A command-line predictor.

The current Nextflow implementation a DSL2 workflow in `main.nf`, 
with the prediction process implemented in
`modules/local/iupred2a/main.nf`. It runs the local `iupred2a/iupred2a.py`
program on one or more FASTA or multi-FASTA files and gives the raw
IUPred2A output files.

## Method

IUPred2A predicts intrinsically disordered protein regions directly from amino
acid sequence. The method estimates whether a sequence segment can form enough
favorable pairwise residue interactions to support a stable globular structure.
Regions with insufficient estimated interaction energy are assigned higher
disorder scores.

This repository also supports ANCHOR2 output through the IUPred2A `-a` option.
ANCHOR2 predicts disordered binding regions, reporting residues that are likely
to be disordered alone but capable of forming interactions with a partner.

Prediction modes:

| Mode | Description |
| --- | --- |
| `long` | Long disordered region prediction. |
| `short` | Short disordered region prediction, such as missing residues in X-ray structures. |
| `glob` | Globular domain prediction. |

For `long` and `short`, the output contains one row per residue with an IUPRED2
score between 0 and 1. Higher values indicate higher predicted disorder, and
scores above 0.5 are commonly treated as disordered. With `--anchor true`, an
additional ANCHOR2 score column is included.

## Requirements

- Nextflow `>=24.04.0`
- Java, as required by Nextflow
- Python 3 when running without a container or Conda profile
- Optional: Docker, Apptainer, or Conda

The bundled IUPred2A Python code does not require external Python packages.

## Quick Start

Run the included example sequence:

```bash
nextflow run doszilab/iupred2a --input 'iupred2a/P53_HUMAN.seq' --mode long
```

Run with ANCHOR2 binding-region prediction:

```bash
nextflow run doszilab/iupred2a --input 'iupred2a/P53_HUMAN.seq' --mode long --anchor true
```

Run all matching FASTA files:

```bash
nextflow run doszilab/iupred2a --input 'data/*.fa' --mode short
```

Run a multi-FASTA file:

```bash
nextflow run doszilab/iupred2a --input 'data/proteins.fasta' --mode long
```

Use Docker:

```bash
nextflow run doszilab/iupred2a --input 'iupred2a/P53_HUMAN.seq' --mode long -profile docker
```

Use Apptainer:

```bash
nextflow run doszilab/iupred2a --input 'iupred2a/P53_HUMAN.seq' --mode long -profile apptainer
```

Use Conda:

```bash
nextflow run doszilab/iupred2a --input 'iupred2a/P53_HUMAN.seq' --mode long -profile conda
```

## Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--input` | required | FASTA or multi-FASTA file path or glob pattern. Each sequence record is processed as one sample. |
| `--mode` | `long` | IUPred2A mode: `long`, `short`, or `glob`. |
| `--anchor` | `false` | Enable ANCHOR2 prediction. Accepts boolean-like values such as `true`, `false`, `1`, `0`, `yes`, or `no`. |
| `--outdir` | `results` | Directory where prediction files are copied. |

## Output

Results are copied to `--outdir`.

Without ANCHOR2:

```text
<sample>.iupred2a.<mode>.txt
```

With ANCHOR2:

```text
<sample>.iupred2a.<mode>.anchor2.txt
```

Single-record FASTA files keep the input file base name as `<sample>`. For
multi-FASTA input, sample names include the input file base name, a zero-padded
record number, and the sanitized first token from the FASTA header:

```text
proteins__000001__first_sequence.iupred2a.long.txt
proteins__000002__second_sequence.iupred2a.long.txt
```

Example output header:

```text
# Prediction type: long
# Prediction output
# POS	RES	IUPRED2
1	M	0.9807
```

With `--anchor true`, the residue table contains:

```text
# POS	RES	IUPRED2	ANCHOR2
```

## Workflow Structure

`main.nf` validates the required input, normalizes `--mode` and `--anchor`, and
creates one Nextflow task per FASTA record. A preprocessing step splits
multi-FASTA inputs into single-sequence FASTA files before prediction. The
`IUPRED2A` process stages each split FASTA and the bundled `iupred2a/`
directory because the Python library loads its `data/` files relative to the
script location.

The command executed by each task is:

```bash
python3 iupred2a/iupred2a.py [-a] <fasta> <mode>
```

## Version Check

This repository was smoke-tested with:

```text
Nextflow 26.04.0
Workflow manifest version 0.1.0
```

The test command used was:

```bash
nextflow run doszilab/iupred2a --input iupred2a/P53_HUMAN.seq --mode long
```

## Citation

If you use IUPred2A results, cite:

Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi. IUPred2A:
context-dependent prediction of protein disorder as a function of redox state
and protein binding. Nucleic Acids Research, 46(W1):W329-W337, 2018.
https://doi.org/10.1093/nar/gky384

For Nextflow, see https://www.nextflow.io/.

## License

The bundled IUPred2A software is distributed with its academic license in
`iupred2a/LICENSE`. Review that license before using or redistributing this
repository, especially for non-academic or commercial use.
