# demo_anarcii

Translates raw nucleotide antibody sequences to amino acids and numbers them with [ANARCII](https://github.com/oxpig/ANARCII) using the IMGT scheme.

## How to

```bash
pip install anarcii pandas numpy
python main.py
```

Reads `input/<file_input>`, numbers every sequence, writes the results table to `output/<file_output>` (paths set in `config.json`).

## File structure

```
main.py                  entry point — reads config.json, processes the input file, writes output CSV
config.json               paths + column name used by main.py
scripts/
  anarcii.py               model setup + process_sequence() (core per-sequence pipeline)
  helpers.py               codon table, translateNT(), read_json()
input/                     source CSV(s)
output/                    generated result CSV(s)
anarcii.ipynb              dev/exploration notebook
.archive/                  earlier notebook versions kept for reference
```

## Config (`config.json`)

| Key | Meaning |
|---|---|
| `dir_input` / `file_input` | Folder + filename of the input CSV (must contain `seq_col`). |
| `dir_output` / `file_output` | Folder + filename the result CSV is written to. |
| `seq_col` | Column in the input CSV holding the raw NT sequences. |
| `limit_rows` | Limit the processing to the first N number of rows set in this line.| 

## Sources

- [ANARCII](https://github.com/oxpig/ANARCII) — the numbering model this pipeline wraps.
- Note: ANARCII can fail on a sequence with `"Model predicted duplicate numbers"`, usually from long runs of `X` (unresolved/`N` bases). `process_sequence(..., remove_x=True)` mitigates this.
