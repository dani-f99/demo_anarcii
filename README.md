# demo_anarcii

Translates raw nucleotide antibody sequences to amino acids and numbers them with [ANARCII](https://github.com/oxpig/ANARCII) using the IMGT scheme, then implants the AnarcII-aligned CDR3 back into the original table.

## How to

```bash
pip install anarcii pandas numpy regex
python main.py
```

Reads the input CSV named in `config.json`, aligns the germline sequence (+ CDR3) with ANARCII, and writes results to `results/<run_name>/` (folder name set by `run_name` in `config.json`).

## File structure

```
main.py                    entry point — reads config.json, runs AssignCDR3, writes results + run report
config.json                 settings + AssignCDR3 column names used by main.py
scripts/
  anarcii.py                 AssignCDR3 class (core alignment/CDR3 pipeline), use_anarchii(), aligned_frame()
  helpers.py                 codon table, translateNT(), get_region(), get_digit(), read_json()
  report.py                  write_report() — timestamped run_report.txt logging
input/                      source CSV(s)
results/<run_name>/         generated per-run outputs (created on run)
  anarcii_sequences.csv       per-position IMGT-aligned residues
  anarcii_output.csv          AnarcII model output + assembled germline/CDR3 strings
  sequences_results.csv       final table (original columns + germline_anarcii + cdr3_anacrii)
  run_report.txt              timestamped log of run steps
dev_notebook.ipynb          dev/exploration notebook
.archived/                  earlier notebook/script versions kept for reference
```

## Config (`config.json`)

| Key | Meaning |
|---|---|
| `settings.run_name` | Name of the run; results are written to `results/<run_name>/`. |
| `settings.dir_input` / `file_input` | Folder + filename of the input CSV. |
| `settings.dir_output` / `file_output` | Reserved output folder/filename settings (current pipeline writes under `results/<run_name>/` instead). |
| `settings.seq_col` | Column in the input CSV holding the raw NT sequences. |
| `settings.limit_rows` | Limit the processing to the first N rows (not currently applied in `main.py`). |
| `AssignCDR3.germline_column` | Column holding the germline NT sequence. |
| `AssignCDR3.sequence_column` | Column holding the sequencing NT sequence. |
| `AssignCDR3.cdr3_column` | Column holding the CDR3 AA sequence to append to the translated germline. |

## Sources

- [ANARCII](https://github.com/oxpig/ANARCII) — the numbering model this pipeline wraps.
- Note: ANARCII can fail on a sequence with `"Model predicted duplicate numbers"`, usually from long runs of `X` (unresolved/`N` bases).
