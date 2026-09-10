# demo_anarcii

Pipeline that takes raw nucleotide (NT) antibody sequences from an ImmuneDB-style export, translates them to amino acids, and numbers them with [ANARCII](https://github.com/oxpig/ANARCII) using the IMGT numbering scheme.

## What it does

1. Reads a CSV of raw NT sequences (with `-` alignment spacers and `N` low-quality/unsequenced bases, as typically produced by NGS reads).
2. Translates each sequence to amino acids in reading frame 1.
3. Strips spacer/gap characters and (optionally) unresolved `X` residues before handing the sequence to ANARCII, since long runs of `X` can cause ANARCII's model to fail with `"Model predicted duplicate numbers"`.
4. Runs ANARCII to assign IMGT positions to each residue.
5. Writes one row per input sequence to an output CSV containing both the raw/original translation and the ANARCII-numbered result, plus ANARCII's own metadata (chain type, score, error, etc).

## Project layout

```
main.py                  entry point — reads config.json, processes the whole input file, writes output CSV
config.json               paths + column name used by main.py
scripts/
  anarcii.py               model setup + process_sequence() (the core per-sequence pipeline)
  helpers.py               codon table, translateNT() (NT -> AA translation), read_json()
input/                     source CSV(s) of raw NT sequences
output/                    generated result CSV(s)
anarcii.ipynb              exploratory/dev notebook used to build and debug the pipeline
.archive/                  earlier notebook versions kept for reference
```

## Setup

```bash
pip install anarcii pandas numpy
```

## Configuration (`config.json`)

```json
{
  "settings": {
    "dir_input": "input",
    "dir_output": "output",
    "file_input": "cleaned_seqs_all_seq_1k.csv",
    "file_output": "1k_test_output.csv",
    "seq_col": "sequence"
  }
}
```

| Key | Meaning |
|---|---|
| `dir_input` / `file_input` | Folder + filename of the input CSV. Must contain the column named by `seq_col`. |
| `dir_output` / `file_output` | Folder + filename the result CSV is written to. |
| `seq_col` | Name of the column in the input CSV holding the raw NT sequence strings. |

## Usage

```bash
python main.py
```

This reads `input/<file_input>`, runs every sequence through `process_sequence(..., anarcii=True, output_type="table")`, and writes the combined results table to `output/<file_output>`.

## Output columns

| Column | Description |
|---|---|
| `nt_og` | Original raw NT sequence (unmodified input). |
| `seq_og` | Plain translation of `nt_og` to amino acids (no ANARCII numbering, spacers/`X` kept as-is). |
| `seq_anarcii` | Amino acid sequence reconstructed from ANARCII's IMGT-numbered output. `NaN` if ANARCII failed to number the sequence. |
| `numbering` | Raw ANARCII output: a list of `((position, insertion_code), amino_acid)` tuples. `None` on failure. |
| `chain_type` | Chain type ANARCII detected (e.g. `H` for heavy, `F` on failure). |
| `score` | ANARCII's confidence score for the numbering. |
| `query_start` / `query_end` | Index range in the translated sequence that ANARCII actually numbered. |
| `error` | ANARCII's error message if numbering failed (e.g. `"Model predicted duplicate numbers"`), otherwise `None`. |
| `scheme` | Numbering scheme used (currently always `"imgt"`). |

## Core functions

### `scripts.helpers.translateNT(seq: str) -> str`
Translates a nucleotide string to amino acids, codon by codon, using the standard genetic code. Unrecognized codons (e.g. containing `N`) translate to `X`; a literal `---` codon translates to `-` (gap).

### `scripts.anarcii.process_sequence(seq, anarcii=True, output_type="str", remove_x=True) -> str | np.ndarray`
The main per-sequence pipeline:
- `anarcii=False` — skip ANARCII entirely, just return the plain translation (`translateNT(seq)`).
- `anarcii=True` — validate that `-` spacer characters occur in complete codon-sized runs, strip them, translate, optionally strip `X` residues (`remove_x`), then run ANARCII.
- `output_type="str"` — return only the ANARCII-numbered amino acid string.
- `output_type="table"` — return a single row (`numpy.ndarray`, in the order of `cols_result` — see the columns table above) combining the original sequence, both translations, and full ANARCII metadata.

**Known limitation:** for a batch of sequences, `.apply(process_sequence, ...)` returns a `Series` where each element is that row-array. Expanding it into a proper DataFrame requires `pd.DataFrame(series.tolist(), columns=cols, index=series.index)` — a plain `pd.DataFrame(series.values)` will *not* expand the arrays into columns.

## Known issues / gotchas

- **`"Model predicted duplicate numbers"`**: ANARCII's own internal sanity check rejects a small fraction of inputs, most often when the translated sequence has long runs of `X` (from `N`-heavy raw reads) confusing its position-prediction model. This is a model limitation, not a bug in this pipeline; `remove_x=True` (the default) mitigates it by stripping `X` residues before calling ANARCII.
- **Spacer validation is total-count-only**: `process_sequence` checks that the *total* number of `-` characters in a sequence is a multiple of 3, not that each individual run of spacers is itself codon-aligned. A sequence with multiple, unevenly-sized spacer runs can pass this check while still having its reading frame corrupted between gap sites.
- **Not all input rows are guaranteed present**: some source CSVs have gaps in their row index (e.g. rows dropped from an earlier filtering step), so positional slicing (`.head(n)`, `.iloc[]`) and label-based indexing (`.loc[]`) can disagree — use `.reset_index()` if you need a clean, dense index.
