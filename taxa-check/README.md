# taxa-check

Check taxon-name consistency between FASTA alignments and Newick trees, clean internal node labels, and batch-rename FASTA headers plus tree leaf labels.

Two entry points are provided:
- `taxa_check.py`: full implementation, runs on Windows / macOS / Linux, supports `check`, `clean`, and `rename`
- `taxa_check.sh`: independent shell implementation for `check` and `clean` only; `rename` is intentionally unsupported

## Usage

Default command is `check`, so old invocations still work:

```bash
# default = check
python taxa_check.py /path/to/data
python taxa_check.py check /path/to/data

# backward-compatible old style
python taxa_check.py --clean /path/to/data
python taxa_check.py --normalize-ext /path/to/data

# explicit subcommands
python taxa_check.py clean /path/to/data
python taxa_check.py rename --rule lower /path/to/data

# shell entry
bash taxa_check.sh check /path/to/data
bash taxa_check.sh clean /path/to/data

# rename is Python-only
python taxa_check.py rename --rule lower /path/to/data
```

Shared flags:
- `--batch`: iterate over first-level subdirectories
- `--logdir DIR`: mismatch log directory name/path, default `taxa_mismatch_logs`
- `--report PATH`: write a TSV report

Shell/Python split:
- `.py`: `check`, `clean`, `rename`
- `.sh`: `check`, `clean`
- `.sh rename`: not supported; use `.py`

`--logdir` behavior:
- relative paths are resolved under the processed target directory
- absolute paths are respected as-is

Supported file extensions:
- FASTA: `.fa` `.fasta` `.faa` `.fna`
- tree: `.nwk` `.tre` `.newick` `.treefile`

## Commands

### `check`

Compare FASTA taxa vs tree leaf labels by stem-matched file pairs.

```bash
python taxa_check.py check --report report.tsv /path/to/data
python taxa_check.py check --clean /path/to/data
python taxa_check.py check --normalize-ext /path/to/data
python taxa_check.py --batch /path/to/parent

bash taxa_check.sh check --report report.tsv /path/to/data
bash taxa_check.sh --clean /path/to/data
bash taxa_check.sh check --normalize-ext /path/to/data
bash taxa_check.sh --batch /path/to/parent
```

Behavior:
- `--normalize-ext` runs before pairing and renames `.newick`, `.treefile`, and `.tre` tree files to `.nwk`
- stem-matched pairs get status `OK` or `MISMATCH`
- unmatched FASTA/tree files get status `NO_TREE` or `NO_FASTA`
- duplicate stems are reported in `notes` as `multiple_fasta=...` / `multiple_tree=...`
- mismatch logdir is created only when a real `MISMATCH` occurs

`--normalize-ext` conflict handling:
- existing `.nwk` files are left unchanged
- if the target `.nwk` path already exists, the source file is skipped with a warning
- `.nex` / `.nexus` files are warned about and skipped because NEXUS is not renamed directly
- reports do not add a normalize column; paths reflect the files after normalization

`check --report` columns:

```text
prefix
status
notes
n_fasta
n_tree
n_tree_only
n_fasta_only
fasta_path
tree_path
log_path
```

### `clean`

Remove internal node labels from tree files only.

```bash
python taxa_check.py clean --report clean.tsv /path/to/data
python taxa_check.py clean --batch /path/to/parent

bash taxa_check.sh clean --report clean.tsv /path/to/data
bash taxa_check.sh clean --batch /path/to/parent
```

Handled label forms include:
- quoted internal labels
- numeric support values
- text node names
- `[&...]` annotations
- IQ-TREE composite support values such as `)100/95:0.05` and `)100/95/0.99:0.05`

Before modification, the script creates `.bak` backups. If `.bak` already exists, a timestamped backup is used.

`clean --report` columns:

```text
prefix
status
n_labels_removed
tree_path
backup_path
```

Statuses:
- `CLEANED`
- `NO_LABELS`
- `NO_TREE`

### `rename`

Rename FASTA headers and tree leaf labels together for each matched prefix.

This command is available in `taxa_check.py` only.

```bash
python taxa_check.py rename --rule lower /path/to/data
python taxa_check.py rename --rule abbrev --parts 3,3 --sep _ --dry-run /path/to/data
python taxa_check.py rename --rule from-file --map mapping.tsv /path/to/data
```

Built-in rules:
- `upper`
- `lower`
- `title` (uppercase first character only, lowercase the rest)
- `abbrev`
- `truncate-dot`
- `join-dot`
- `replace`
- `strip-version`
- `from-file`

Rule-specific options:
- `abbrev`: `--parts`, `--sep` (`--sep ""` produces names like `HomSap`)
- `truncate-dot`: `--keep-parts`
- `replace`: `--from`, `--to`
- `from-file`: `--map`, optional `--strict-mapping`

Rename behavior:
- writes `<prefix>_rename_mapping.tsv` next to the source files
- `--dry-run` writes mapping only, without modifying FASTA/tree files
- apply mode creates `.bak` backups for both FASTA and tree files
- collisions are detected before writing and abort that prefix
- `--strict-mapping` failures are reported as `ERROR`

`rename --report` columns:

```text
prefix
status
n_renamed
n_unchanged
n_collisions
mapping_path
```

Statuses:
- `APPLIED`
- `DRY_RUN`
- `COLLISION`
- `ERROR`
- `NO_TREE`
- `NO_FASTA`
- `NO_CHANGE`

## Batch Mode

With `--batch`, the tool scans only the first level of subdirectories under the given parent directory:

```text
parent/
  gene_A/
  gene_B/
  ...
```

Rules:
- subdirectories are processed in lexical order
- a child directory whose name equals `--logdir` is skipped
- report rows are merged into one TSV
- `prefix` is written as `subdir/prefix` to avoid collisions

## Exit Codes

- `check`: number of mismatched pairs, capped at 255
- `clean`: number of `NO_TREE` rows, capped at 255
- `rename`: collision/strict-mapping failures accumulate in non-zero exit status, capped at 255

## Notes

- FASTA taxa are parsed as the first whitespace-delimited token after `>`
- quoted tree labels preserve spaces
- doubled quotes inside quoted Newick labels are decoded correctly, e.g. `'O''Brien'`
- UTF-8 BOM input is handled automatically
- `taxa_check.sh` is useful on shell-only systems where Python may be unavailable, but it supports only `check` and `clean`
