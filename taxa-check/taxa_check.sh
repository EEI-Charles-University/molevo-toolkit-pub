#!/usr/bin/env bash
# ============================================================================
#  taxa_check.sh - independent shell implementation for check/clean only
#                  rename is intentionally unsupported here; use taxa_check.py
# ============================================================================

set -euo pipefail
shopt -s nullglob

FASTA_EXTS=".fa .fasta .faa .fna"
TREE_EXTS=".nwk .tre .newick .treefile"
NORMALIZE_TREE_EXTS=".newick .treefile .tre"
NEXUS_EXTS=".nex .nexus"

COMMAND="check"
BATCH=false
CHECK_CLEAN=false
CHECK_NORMALIZE_EXT=false
LOGDIR="taxa_mismatch_logs"
REPORT_PATH=""
TARGET_DIR="."
REPORT_TMP=""

declare -a _TMPFILES=()
declare -a PREFIXES=()
declare -A FA_BY=()
declare -A TR_BY=()
declare -A NOTE_BY=()

_cleanup_tmpfiles() {
  if (( ${#_TMPFILES[@]} > 0 )); then
    rm -f -- "${_TMPFILES[@]}" 2>/dev/null || true
  fi
}

trap _cleanup_tmpfiles EXIT

_print_help() {
  cat <<'EOF'
Batch check taxon name consistency between FASTA and tree files
批量检查 FASTA 与树文件之间的物种名一致性（Shell 基础版）

Usage:
  taxa_check.sh [check] [OPTIONS] [DIRECTORY]
  taxa_check.sh clean   [OPTIONS] [DIRECTORY]

Supported in .sh:
  check    consistency check, optional --clean, optional --report
  clean    internal-label cleaning only, optional --report

Not supported in .sh:
  rename   please use taxa_check.py rename ...

Shared options:
  --batch          iterate first-level subdirectories
  --logdir DIR     mismatch log directory (default: taxa_mismatch_logs)
  --report PATH    write TSV report
  -h, --help       show this help

check-only option:
  --clean          clean internal node labels before check
  --normalize-ext  rename .newick/.treefile/.tre tree files to .nwk before check

Examples:
  bash taxa_check.sh /path/to/data
  bash taxa_check.sh check --report report.tsv /path/to/data
  bash taxa_check.sh --clean /path/to/data
  bash taxa_check.sh clean --batch /path/to/parent

Note:
  rename is intentionally unavailable in shell mode.
  Use: python taxa_check.py rename ...
EOF
}

_die() {
  echo "$1" >&2
  exit "${2:-2}"
}

_append_note() {
  local key="$1"
  local note="$2"
  if [[ -n "${NOTE_BY[$key]:-}" ]]; then
    NOTE_BY["$key"]="${NOTE_BY[$key]};${note}"
  else
    NOTE_BY["$key"]="$note"
  fi
}

_format_prefix() {
  local prefix_dir="$1"
  local prefix="$2"
  if [[ -n "$prefix_dir" ]]; then
    printf '%s/%s' "$prefix_dir" "$prefix"
  else
    printf '%s' "$prefix"
  fi
}

_is_abs_path() {
  local path="$1"
  [[ "$path" == /* || "$path" =~ ^[A-Za-z]:[\\/].* ]]
}

_resolve_logdir() {
  local directory="$1"
  local logdir="$2"
  if _is_abs_path "$logdir"; then
    printf '%s' "$logdir"
  else
    printf '%s/%s' "$directory" "$logdir"
  fi
}

_safe_backup() {
  local filepath="$1"
  local backup="${filepath}.bak"
  if [[ -e "$backup" ]]; then
    local ts
    ts="$(date +%Y%m%d_%H%M%S)"
    backup="${filepath}.bak_${ts}"
  fi
  cp "$filepath" "$backup"
  printf '%s' "$backup"
}

_track_tmp() {
  local path="$1"
  _TMPFILES+=("$path")
}

_untrack_tmp() {
  local path="$1"
  local -a kept=()
  local item
  for item in "${_TMPFILES[@]}"; do
    [[ "$item" == "$path" ]] || kept+=("$item")
  done
  _TMPFILES=("${kept[@]}")
}

_mktemp_tracked() {
  local tmp
  tmp="$(mktemp)"
  _track_tmp "$tmp"
  printf '%s' "$tmp"
}

_rm_tmp() {
  local path="$1"
  rm -f -- "$path"
  _untrack_tmp "$path"
}

_report_init() {
  local header="$1"
  if [[ -z "$REPORT_PATH" ]]; then
    REPORT_TMP=""
    return
  fi
  REPORT_TMP="$(_mktemp_tracked)"
  printf '%s\n' "$header" > "$REPORT_TMP"
}

_report_append() {
  if [[ -z "$REPORT_TMP" ]]; then
    return
  fi
  printf '%s\n' "$1" >> "$REPORT_TMP"
}

_report_finalize() {
  if [[ -z "$REPORT_TMP" ]]; then
    return
  fi
  cp "$REPORT_TMP" "$REPORT_PATH"
  _rm_tmp "$REPORT_TMP"
  REPORT_TMP=""
}

_fasta_taxa() {
  awk '
    NR == 1 { sub(/^\357\273\277/, "", $0) }
    /^>/ {
      line = substr($0, 2)
      sub(/^[ \t]+/, "", line)
      if (line != "") {
        split(line, parts, /[ \t]/)
        if (parts[1] != "") print parts[1]
      }
    }
  ' "$1" | LC_ALL=C sort -u
}

_tree_taxa() {
  tr -d '\r\n' < "$1" | awk '
    function read_quoted(quote,    ch, label) {
      label = ""
      i++
      while (i <= n) {
        ch = substr($0, i, 1)
        if (ch == quote) {
          if (i < n && substr($0, i + 1, 1) == quote) {
            label = label quote
            i += 2
            continue
          }
          i++
          return label
        }
        label = label ch
        i++
      }
      return label
    }
    BEGIN { expect_leaf = 1 }
    {
      sub(/^\357\273\277/, "", $0)
      n = length($0)
      i = 1
      while (i <= n) {
        ch = substr($0, i, 1)
        if (ch == "(") {
          expect_leaf = 1
          i++
        } else if (ch == ")") {
          expect_leaf = 0
          i++
        } else if (ch == ",") {
          expect_leaf = 1
          i++
        } else if (ch == ";") {
          i++
        } else if (ch == "[") {
          depth = 1
          i++
          while (i <= n && depth > 0) {
            c2 = substr($0, i, 1)
            if (c2 == "[") depth++
            else if (c2 == "]") depth--
            i++
          }
        } else if (ch == ":") {
          i++
          while (i <= n) {
            c2 = substr($0, i, 1)
            if (index("0123456789.eE+-", c2) == 0) break
            i++
          }
        } else if (ch == "'"'"'" || ch == "\"") {
          label = read_quoted(ch)
          if (expect_leaf && label != "") print label
        } else if (ch == " " || ch == "\t") {
          i++
        } else {
          start = i
          while (i <= n) {
            c2 = substr($0, i, 1)
            if (index(":,();[] \t", c2) > 0) break
            i++
          }
          label = substr($0, start, i - start)
          gsub(/^[ \t]+|[ \t]+$/, "", label)
          if (expect_leaf && label != "") print label
        }
      }
    }
  ' | LC_ALL=C sort -u
}

_clean_tree_to_file() {
  local input="$1"
  local output="$2"
  local count_file="$3"

  tr -d '\r\n' < "$input" | awk -v count_file="$count_file" '
    function read_annotation(    ch, text, depth) {
      text = "["
      depth = 1
      i++
      while (i <= n && depth > 0) {
        ch = substr($0, i, 1)
        text = text ch
        if (ch == "[") depth++
        else if (ch == "]") depth--
        i++
      }
      return text
    }
    function read_quoted_raw(quote,    ch, raw) {
      raw = quote
      i++
      while (i <= n) {
        ch = substr($0, i, 1)
        raw = raw ch
        if (ch == quote) {
          if (i < n && substr($0, i + 1, 1) == quote) {
            raw = raw quote
            i += 2
            continue
          }
          i++
          return raw
        }
        i++
      }
      return raw
    }
    BEGIN { removed = 0; out = ""; expect_leaf = 1 }
    {
      sub(/^\357\273\277/, "", $0)
      n = length($0)
      i = 1
      while (i <= n) {
        ch = substr($0, i, 1)
        if (ch == "(") {
          out = out ch
          expect_leaf = 1
          i++
        } else if (ch == ")") {
          out = out ch
          expect_leaf = 0
          i++
        } else if (ch == ",") {
          out = out ch
          expect_leaf = 1
          i++
        } else if (ch == ";") {
          out = out ch
          i++
        } else if (ch == "[") {
          ann = read_annotation()
          if (expect_leaf) out = out ann
          else removed++
        } else if (ch == ":") {
          start = i
          i++
          while (i <= n) {
            c2 = substr($0, i, 1)
            if (index("0123456789.eE+-", c2) == 0) break
            i++
          }
          out = out substr($0, start, i - start)
        } else if (ch == "'"'"'" || ch == "\"") {
          raw = read_quoted_raw(ch)
          if (expect_leaf) out = out raw
          else removed++
        } else if (ch == " " || ch == "\t") {
          out = out ch
          i++
        } else {
          start = i
          while (i <= n) {
            c2 = substr($0, i, 1)
            if (index(":,();[] \t", c2) > 0) break
            i++
          }
          token = substr($0, start, i - start)
          if (expect_leaf) out = out token
          else if (token != "") removed++
        }
      }
    }
    END {
      print out
      print removed > count_file
    }
  ' > "$output"
}

CLEAN_LAST_REMOVED=0
CLEAN_LAST_BACKUP=""
_clean_internal_labels() {
  local filepath="$1"
  local tmp_out tmp_count removed backup
  tmp_out="$(_mktemp_tracked)"
  tmp_count="$(_mktemp_tracked)"

  _clean_tree_to_file "$filepath" "$tmp_out" "$tmp_count"
  removed="$(tr -d '[:space:]' < "$tmp_count")"
  _rm_tmp "$tmp_count"
  if [[ -z "$removed" ]]; then
    removed=0
  fi

  if (( removed > 0 )); then
    backup="$(_safe_backup "$filepath")"
    cp "$tmp_out" "$filepath"
    CLEAN_LAST_BACKUP="$backup"
    CLEAN_LAST_REMOVED="$removed"
  else
    CLEAN_LAST_BACKUP=""
    CLEAN_LAST_REMOVED=0
  fi

  _rm_tmp "$tmp_out"
}

_normalize_tree_extensions() {
  local directory="$1"
  local file base stem ext target_base target_path
  local renamed=0 candidates=0

  echo "🔧 统一树文件扩展名为 .nwk... / Normalizing tree file extensions to .nwk..."
  for file in "$directory"/*; do
    [[ -f "$file" ]] || continue
    base="$(basename "$file")"
    stem="${base%.*}"
    ext=".${base##*.}"
    ext="${ext,,}"

    if [[ " $NEXUS_EXTS " == *" $ext "* ]]; then
      printf '⚠️  %s: NEXUS 格式不支持直接重命名，跳过\n' "$base"
      printf '⚠️  %s: NEXUS format not supported for renaming, skipped\n' "$base"
      continue
    fi

    if [[ " $NORMALIZE_TREE_EXTS " != *" $ext "* ]]; then
      continue
    fi

    candidates=$((candidates + 1))
    target_base="${stem}.nwk"
    target_path="${directory}/${target_base}"

    if [[ -e "$target_path" ]]; then
      printf '⚠️  %s: 无法重命名 %s → %s，目标文件已存在，跳过\n' "$stem" "$base" "$target_base"
      printf '⚠️  %s: Cannot rename %s → %s, target already exists, skipped\n' "$stem" "$base" "$target_base"
      continue
    fi

    mv -- "$file" "$target_path"
    printf '   renamed 已重命名: %s → %s\n' "$base" "$target_base"
    renamed=$((renamed + 1))
  done

  if (( candidates == 0 )); then
    echo "   所有树文件已是 .nwk / All tree files already use .nwk"
  else
    printf '   共重命名 %d 个文件 / Renamed %d file(s)\n' "$renamed" "$renamed"
  fi
  echo
}

_collect_files() {
  local directory="$1"
  local file base stem ext used ignored others
  local -A fa_dups=()
  local -A tr_dups=()

  FA_BY=()
  TR_BY=()
  NOTE_BY=()
  PREFIXES=()

  for file in "$directory"/*; do
    [[ -f "$file" ]] || continue
    base="$(basename "$file")"
    stem="${base%.*}"
    ext=".${base##*.}"
    ext="${ext,,}"

    if [[ " $FASTA_EXTS " == *" $ext "* ]]; then
      if [[ -n "${fa_dups[$stem]:-}" ]]; then
        fa_dups["$stem"]="${fa_dups[$stem]}|$base"
      else
        fa_dups["$stem"]="$base"
      fi
      [[ -n "${FA_BY[$stem]:-}" ]] || FA_BY["$stem"]="$file"
    elif [[ " $TREE_EXTS " == *" $ext "* ]]; then
      if [[ -n "${tr_dups[$stem]:-}" ]]; then
        tr_dups["$stem"]="${tr_dups[$stem]}|$base"
      else
        tr_dups["$stem"]="$base"
      fi
      [[ -n "${TR_BY[$stem]:-}" ]] || TR_BY["$stem"]="$file"
    fi
  done

  local warned=false
  while IFS= read -r stem; do
    [[ -n "$stem" ]] || continue
    IFS='|' read -r -a names <<< "${fa_dups[$stem]}"
    if (( ${#names[@]} > 1 )); then
      _append_note "$stem" "multiple_fasta=${fa_dups[$stem]}"
      used="$(basename "${FA_BY[$stem]}")"
      others=""
      for base in "${names[@]}"; do
        [[ "$base" == "$used" ]] && continue
        if [[ -n "$others" ]]; then
          others="${others}, ${base}"
        else
          others="$base"
        fi
      done
      printf '⚠️  %-25s MULTIPLE_FASTA: using 使用 %s, ignoring 忽略 %s\n' \
        "$stem" "$used" "$others"
      warned=true
    fi
  done < <(printf '%s\n' "${!fa_dups[@]}" | LC_ALL=C sort)

  while IFS= read -r stem; do
    [[ -n "$stem" ]] || continue
    IFS='|' read -r -a names <<< "${tr_dups[$stem]}"
    if (( ${#names[@]} > 1 )); then
      _append_note "$stem" "multiple_tree=${tr_dups[$stem]}"
      used="$(basename "${TR_BY[$stem]}")"
      others=""
      for base in "${names[@]}"; do
        [[ "$base" == "$used" ]] && continue
        if [[ -n "$others" ]]; then
          others="${others}, ${base}"
        else
          others="$base"
        fi
      done
      printf '⚠️  %-25s MULTIPLE_TREE:  using 使用 %s, ignoring 忽略 %s\n' \
        "$stem" "$used" "$others"
      warned=true
    fi
  done < <(printf '%s\n' "${!tr_dups[@]}" | LC_ALL=C sort)

  if $warned; then
    echo
  fi

  mapfile -t PREFIXES < <(
    {
      printf '%s\n' "${!FA_BY[@]}"
      printf '%s\n' "${!TR_BY[@]}"
    } | awk 'NF && !seen[$0]++ { print $0 }' | LC_ALL=C sort
  )
}

CHECK_BAD=0
_run_check_directory() {
  local directory="$1"
  local logdir="$2"
  local prefix_dir="$3"
  local prefix display_prefix notes fasta_path tree_path log_path
  local matched=0 ok=0 bad=0
  local fa_taxa tree_taxa
  local t1 t2 n_fasta n_tree n_tree_only n_fasta_only

  if $CHECK_NORMALIZE_EXT; then
    _normalize_tree_extensions "$directory"
  fi

  _collect_files "$directory"

  if $CHECK_CLEAN; then
    echo "🔧 清理树文件中的内部节点标签... / Cleaning internal node labels..."
    local cleaned=0
    while IFS= read -r prefix; do
      [[ -n "$prefix" ]] || continue
      _clean_internal_labels "${TR_BY[$prefix]}"
      if (( CLEAN_LAST_REMOVED > 0 )); then
        printf '   cleaned 已清理: %s (backup 备份: %s)\n' \
          "${TR_BY[$prefix]}" "$CLEAN_LAST_BACKUP"
        cleaned=$((cleaned + 1))
      fi
    done < <(printf '%s\n' "${!TR_BY[@]}" | LC_ALL=C sort)
    if (( cleaned == 0 )); then
      echo "   未发现需要清理的内部节点标签 / No internal node labels found"
    else
      printf '   共清理 %d 个树文件 / Cleaned %d tree file(s)\n' "$cleaned" "$cleaned"
    fi
    echo
  fi

  if (( ${#PREFIXES[@]} == 0 )); then
    printf '⚠️  目录为空或无可识别文件 / No FASTA or tree files found: %s\n' "$directory"
    echo "----"
    echo "Paired 已配对: 0   OK 通过: 0   MISMATCH 不一致: 0"
    CHECK_BAD=0
    return
  fi

  for prefix in "${PREFIXES[@]}"; do
    display_prefix="$(_format_prefix "$prefix_dir" "$prefix")"
    notes="${NOTE_BY[$prefix]:-}"
    fasta_path="${FA_BY[$prefix]:-}"
    tree_path="${TR_BY[$prefix]:-}"

    if [[ -n "$fasta_path" && -n "$tree_path" ]]; then
      matched=$((matched + 1))
      t1="$(_mktemp_tracked)"
      t2="$(_mktemp_tracked)"
      _fasta_taxa "$fasta_path" > "$t1"
      _tree_taxa "$tree_path" > "$t2"

      n_fasta="$(wc -l < "$t1" | tr -d ' ')"
      n_tree="$(wc -l < "$t2" | tr -d ' ')"
      n_tree_only="$(LC_ALL=C comm -23 "$t2" "$t1" | wc -l | tr -d ' ')"
      n_fasta_only="$(LC_ALL=C comm -13 "$t2" "$t1" | wc -l | tr -d ' ')"

      if [[ "$n_tree_only" == "0" && "$n_fasta_only" == "0" ]]; then
        printf '✅ %-25s OK (n=%s)\n' "$prefix" "$n_tree"
        ok=$((ok + 1))
        log_path=""
        _report_append "${display_prefix}	OK	${notes}	${n_fasta}	${n_tree}	0	0	${fasta_path}	${tree_path}	"
      else
        mkdir -p "$logdir"
        log_path="${logdir}/${prefix}.mismatch.txt"
        {
          printf 'PREFIX: %s\n' "$prefix"
          printf 'FASTA:  %s\n' "$fasta_path"
          printf 'TREE:   %s\n' "$tree_path"
          printf 'Counts: TREE=%s FASTA=%s\n\n' "$n_tree" "$n_fasta"
          printf 'In TREE but NOT in FASTA / 在树中但不在FASTA中 (%s):\n' "$n_tree_only"
          LC_ALL=C comm -23 "$t2" "$t1" | sed 's/^/  /'
          printf '\nIn FASTA but NOT in TREE / 在FASTA中但不在树中 (%s):\n' "$n_fasta_only"
          LC_ALL=C comm -13 "$t2" "$t1" | sed 's/^/  /'
          cat <<'EOF'

Hint 提示: also check for case / underscore / dot / dash differences if counts look close.
如果数量接近，请检查大小写、下划线、点、横线等格式差异。

Note 注意: FASTA taxa = first token after '>' (space-delimited); tree may use quoted names with spaces.
FASTA 物种名 = '>' 后第一个空白前的部分；树中引号包裹的名称会保留空格。
EOF
        } > "$log_path"
        printf '❌ %-25s MISMATCH 不一致 (tree_only=%s fasta_only=%s) -> %s\n' \
          "$prefix" "$n_tree_only" "$n_fasta_only" "$log_path"
        bad=$((bad + 1))
        _report_append "${display_prefix}	MISMATCH	${notes}	${n_fasta}	${n_tree}	${n_tree_only}	${n_fasta_only}	${fasta_path}	${tree_path}	${log_path}"
      fi

      _rm_tmp "$t1"
      _rm_tmp "$t2"
      continue
    fi

    if [[ -n "$fasta_path" ]]; then
      t1="$(_mktemp_tracked)"
      _fasta_taxa "$fasta_path" > "$t1"
      n_fasta="$(wc -l < "$t1" | tr -d ' ')"
      _rm_tmp "$t1"
      printf '⚠️  %-25s NO_TREE 缺树文件  (have FASTA 有FASTA: %s)\n' "$prefix" "$fasta_path"
      _report_append "${display_prefix}	NO_TREE	${notes}	${n_fasta}	0	0	${n_fasta}	${fasta_path}		"
    else
      t2="$(_mktemp_tracked)"
      _tree_taxa "$tree_path" > "$t2"
      n_tree="$(wc -l < "$t2" | tr -d ' ')"
      _rm_tmp "$t2"
      printf '⚠️  %-25s NO_FASTA 缺FASTA (have TREE 有树文件: %s)\n' "$prefix" "$tree_path"
      _report_append "${display_prefix}	NO_FASTA	${notes}	0	${n_tree}	${n_tree}	0		${tree_path}	"
    fi
  done

  echo "----"
  printf 'Paired 已配对: %d   OK 通过: %d   MISMATCH 不一致: %d\n' "$matched" "$ok" "$bad"
  if (( bad > 0 )); then
    printf 'Mismatch logs 不一致日志: %s/\n' "$logdir"
  fi
  CHECK_BAD="$bad"
}

CLEAN_FAILS=0
_run_clean_directory() {
  local directory="$1"
  local prefix_dir="$2"
  local prefix display_prefix tree_path fasta_path
  local cleaned=0 no_tree=0

  _collect_files "$directory"

  echo "🔧 清理树文件中的内部节点标签... / Cleaning internal node labels..."

  if (( ${#PREFIXES[@]} == 0 )); then
    echo "   未发现树文件 / No tree files found"
    CLEAN_FAILS=0
    return
  fi

  for prefix in "${PREFIXES[@]}"; do
    display_prefix="$(_format_prefix "$prefix_dir" "$prefix")"
    tree_path="${TR_BY[$prefix]:-}"
    fasta_path="${FA_BY[$prefix]:-}"

    if [[ -z "$tree_path" ]]; then
      printf '⚠️  %-25s NO_TREE 缺树文件  (have FASTA 有FASTA: %s)\n' "$prefix" "$fasta_path"
      _report_append "${display_prefix}	NO_TREE	0		"
      no_tree=$((no_tree + 1))
      continue
    fi

    _clean_internal_labels "$tree_path"
    if (( CLEAN_LAST_REMOVED > 0 )); then
      printf '   cleaned 已清理: %s (backup 备份: %s)\n' "$tree_path" "$CLEAN_LAST_BACKUP"
      _report_append "${display_prefix}	CLEANED	${CLEAN_LAST_REMOVED}	${tree_path}	${CLEAN_LAST_BACKUP}"
      cleaned=$((cleaned + 1))
    else
      _report_append "${display_prefix}	NO_LABELS	0	${tree_path}	"
    fi
  done

  if (( cleaned == 0 )); then
    echo "   未发现需要清理的内部节点标签 / No internal node labels found"
  else
    printf '   共清理 %d 个树文件 / Cleaned %d tree file(s)\n' "$cleaned" "$cleaned"
  fi

  CLEAN_FAILS="$no_tree"
}

_run_check() {
  local total_bad=0
  local subdir full_path prefix_dir resolved_logdir

  _report_init $'prefix\tstatus\tnotes\tn_fasta\tn_tree\tn_tree_only\tn_fasta_only\tfasta_path\ttree_path\tlog_path'

  if $BATCH; then
    local found=false
    while IFS= read -r subdir; do
      [[ -n "$subdir" ]] || continue
      found=true
      full_path="${TARGET_DIR}/${subdir}"
      prefix_dir="$subdir"
      resolved_logdir="$(_resolve_logdir "$full_path" "$LOGDIR")"
      printf '===== %s =====\n' "$subdir"
      _run_check_directory "$full_path" "$resolved_logdir" "$prefix_dir"
      total_bad=$((total_bad + CHECK_BAD))
      echo
    done < <(
      for full_path in "$TARGET_DIR"/*; do
        [[ -d "$full_path" ]] || continue
        subdir="$(basename "$full_path")"
        [[ "$subdir" == "$LOGDIR" ]] && continue
        printf '%s\n' "$subdir"
      done | LC_ALL=C sort
    )

    if ! $found; then
      _report_finalize
      _die "未找到子目录 / No subdirectories found: $TARGET_DIR" 1
    fi
  else
    resolved_logdir="$(_resolve_logdir "$TARGET_DIR" "$LOGDIR")"
    _run_check_directory "$TARGET_DIR" "$resolved_logdir" ""
    total_bad="$CHECK_BAD"
  fi

  _report_finalize
  if (( total_bad > 255 )); then
    total_bad=255
  fi
  exit "$total_bad"
}

_run_clean() {
  local total_fails=0
  local subdir full_path

  _report_init $'prefix\tstatus\tn_labels_removed\ttree_path\tbackup_path'

  if $BATCH; then
    local found=false
    while IFS= read -r subdir; do
      [[ -n "$subdir" ]] || continue
      found=true
      full_path="${TARGET_DIR}/${subdir}"
      printf '===== %s =====\n' "$subdir"
      _run_clean_directory "$full_path" "$subdir"
      total_fails=$((total_fails + CLEAN_FAILS))
      echo
    done < <(
      for full_path in "$TARGET_DIR"/*; do
        [[ -d "$full_path" ]] || continue
        subdir="$(basename "$full_path")"
        [[ "$subdir" == "$LOGDIR" ]] && continue
        printf '%s\n' "$subdir"
      done | LC_ALL=C sort
    )

    if ! $found; then
      _report_finalize
      _die "未找到子目录 / No subdirectories found: $TARGET_DIR" 1
    fi
  else
    _run_clean_directory "$TARGET_DIR" ""
    total_fails="$CLEAN_FAILS"
  fi

  _report_finalize
  if (( total_fails > 255 )); then
    total_fails=255
  fi
  exit "$total_fails"
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    check|clean)
      COMMAND="$1"
      shift
      ;;
    rename)
      _die "rename is not supported in taxa_check.sh; please use: python taxa_check.py rename ..." 2
      ;;
    -h|--help)
      _print_help
      exit 0
      ;;
  esac
fi

declare -a POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean)
      CHECK_CLEAN=true
      shift
      ;;
    --normalize-ext)
      CHECK_NORMALIZE_EXT=true
      shift
      ;;
    --batch)
      BATCH=true
      shift
      ;;
    --logdir)
      [[ $# -ge 2 ]] || _die "--logdir requires a value" 2
      LOGDIR="$2"
      shift 2
      ;;
    --report)
      [[ $# -ge 2 ]] || _die "--report requires a value" 2
      REPORT_PATH="$2"
      shift 2
      ;;
    rename)
      _die "rename is not supported in taxa_check.sh; please use: python taxa_check.py rename ..." 2
      ;;
    -h|--help)
      _print_help
      exit 0
      ;;
    -*)
      _die "Unknown option / 未知选项: $1" 2
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

if (( ${#POSITIONAL[@]} > 0 )); then
  TARGET_DIR="${POSITIONAL[0]}"
fi

if [[ "$COMMAND" == "clean" && "$CHECK_CLEAN" == true ]]; then
  _die "--clean is only meaningful for the check command" 2
fi

if [[ "$COMMAND" == "clean" && "$CHECK_NORMALIZE_EXT" == true ]]; then
  _die "--normalize-ext is only meaningful for the check command" 2
fi

if [[ "$COMMAND" == "check" ]]; then
  _run_check
elif [[ "$COMMAND" == "clean" ]]; then
  _run_clean
else
  _die "Unsupported command: $COMMAND" 2
fi
