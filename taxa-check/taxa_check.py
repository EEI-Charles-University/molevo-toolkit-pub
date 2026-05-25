#!/usr/bin/env python3
# ============================================================================
#  taxa_check.py — 批量检查 FASTA 与树文件之间的物种名一致性，并提供清理/重命名功能
#                  Batch check FASTA/tree taxon consistency, with clean/rename.
# ============================================================================

import argparse
import csv
import os
import re
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path


FASTA_EXTS = {".fa", ".fasta", ".faa", ".fna"}
TREE_EXTS = {".nwk", ".tre", ".newick", ".treefile"}
NORMALIZE_TREE_EXTS = {".newick", ".treefile", ".tre"}
NEXUS_EXTS = {".nex", ".nexus"}
CHECK_REPORT_FIELDS = [
    "prefix",
    "status",
    "notes",
    "n_fasta",
    "n_tree",
    "n_tree_only",
    "n_fasta_only",
    "fasta_path",
    "tree_path",
    "log_path",
]
CLEAN_REPORT_FIELDS = [
    "prefix",
    "status",
    "n_labels_removed",
    "tree_path",
    "backup_path",
]
RENAME_REPORT_FIELDS = [
    "prefix",
    "status",
    "n_renamed",
    "n_unchanged",
    "n_collisions",
    "mapping_path",
]
UNQUOTED_TREE_LABEL_RE = re.compile(r"^[^:\s,();\[\]'\"`]+$")


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")


@dataclass
class FileCollection:
    fasta_map: dict = field(default_factory=dict)
    tree_map: dict = field(default_factory=dict)
    note_map: dict = field(default_factory=lambda: defaultdict(list))
    warning_messages: list = field(default_factory=list)


def _read_text(filepath):
    data = Path(filepath).read_bytes()
    text = data.decode("utf-8", errors="replace")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return text.lstrip("\ufeff")


def _write_text(filepath, text):
    with open(filepath, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _safe_backup(filepath):
    backup = filepath + ".bak"
    if os.path.exists(backup):
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = f"{filepath}.bak_{ts}"
    shutil.copy2(filepath, backup)
    return backup


def _format_prefix(prefix_dir, prefix):
    return f"{prefix_dir}/{prefix}" if prefix_dir else prefix


def _write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def extract_fasta_taxa(filepath):
    taxa = set()
    text = _read_text(filepath)
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0] if line[1:].strip() else ""
            if name:
                taxa.add(name)
    return sorted(taxa)


def _consume_quoted_label(text, start, quote_char):
    i = start + 1
    n = len(text)
    chars = []

    while i < n:
        ch = text[i]
        if ch == quote_char:
            if i + 1 < n and text[i + 1] == quote_char:
                chars.append(quote_char)
                i += 2
                continue
            return "".join(chars), i + 1
        chars.append(ch)
        i += 1

    return "".join(chars), i


def extract_tree_taxa(filepath):
    text = _read_text(filepath).replace("\n", "").strip()
    if not text:
        return []

    taxa = set()
    i = 0
    n = len(text)
    expect_leaf = True

    while i < n:
        ch = text[i]

        if ch == "(":
            expect_leaf = True
            i += 1
        elif ch == ")":
            expect_leaf = False
            i += 1
        elif ch == ",":
            expect_leaf = True
            i += 1
        elif ch == ";":
            i += 1
        elif ch == "[":
            depth = 1
            i += 1
            while i < n and depth > 0:
                if text[i] == "[":
                    depth += 1
                elif text[i] == "]":
                    depth -= 1
                i += 1
        elif ch == ":":
            i += 1
            while i < n and text[i] in "0123456789.eE+-":
                i += 1
        elif ch == "'":
            label, i = _consume_quoted_label(text, i, "'")
            if expect_leaf and label:
                taxa.add(label)
        elif ch == '"':
            label, i = _consume_quoted_label(text, i, '"')
            if expect_leaf and label:
                taxa.add(label)
        elif ch in " \t":
            i += 1
        else:
            start = i
            while i < n and text[i] not in ":,();[] \t":
                i += 1
            label = text[start:i].strip()
            if expect_leaf and label:
                taxa.add(label)

    return sorted(taxa)


_CLEAN_PATTERNS = [
    (re.compile(r"\)'[^']*':"), "):"),
    (re.compile(r'\)"[^"]*":'), "):"),
    (re.compile(r"\)'[^']*'([),;\n])"), r")\1"),
    (re.compile(r'\)"[^"]*"([),;\n])'), r")\1"),
    (re.compile(r"\)'[^']*$", re.MULTILINE), ")"),
    (re.compile(r'\)"[^"]*"$', re.MULTILINE), ")"),
    (re.compile(r"\)\[&[^\]]*\]"), ")"),
    (
        re.compile(r"\)(\d+(?:\.\d+)?)(?:/\d+(?:\.\d+)?){1,2}([:),;])"),
        r")\2",
    ),
    (
        re.compile(
            r"\)(\d+(?:\.\d+)?)(?:/\d+(?:\.\d+)?){1,2}$",
            re.MULTILINE,
        ),
        ")",
    ),
    (re.compile(r"\)([A-Za-z_][A-Za-z0-9_.-]*)([:),;])"), r")\2"),
    (re.compile(r"\)([A-Za-z_][A-Za-z0-9_.-]*)$", re.MULTILINE), ")"),
    (re.compile(r"\)(\d+(?:\.\d+)?)([:),;])"), r")\2"),
    (re.compile(r"\)(\d+(?:\.\d+)?)$", re.MULTILINE), ")"),
]


def has_internal_labels(filepath):
    text = _read_text(filepath)
    return any(pattern.search(text) for pattern, _ in _CLEAN_PATTERNS)


def clean_tree_text(text):
    total = 0
    new_text = text
    for pattern, replacement in _CLEAN_PATTERNS:
        new_text, count = pattern.subn(replacement, new_text)
        total += count
    return new_text, total


def clean_internal_labels(filepath):
    text = _read_text(filepath)
    cleaned_text, removed = clean_tree_text(text)
    if removed == 0:
        return "", 0

    backup = _safe_backup(filepath)
    _write_text(filepath, cleaned_text)
    return backup, removed


def normalize_tree_extensions(directory):
    print("🔧 统一树文件扩展名为 .nwk... / Normalizing tree file extensions to .nwk...")
    renamed = 0
    candidates = 0

    for entry in sorted(os.listdir(directory)):
        full_path = os.path.join(directory, entry)
        if not os.path.isfile(full_path):
            continue

        stem, ext = os.path.splitext(entry)
        ext = ext.lower()

        if ext in NEXUS_EXTS:
            print(f"⚠️  {entry}: NEXUS 格式不支持直接重命名，跳过")
            print(f"⚠️  {entry}: NEXUS format not supported for renaming, skipped")
            continue

        if ext not in NORMALIZE_TREE_EXTS:
            continue

        candidates += 1
        target_name = f"{stem}.nwk"
        target_path = os.path.join(directory, target_name)

        if os.path.exists(target_path):
            print(f"⚠️  {stem}: 无法重命名 {entry} → {target_name}，目标文件已存在，跳过")
            print(f"⚠️  {stem}: Cannot rename {entry} → {target_name}, target already exists, skipped")
            continue

        os.rename(full_path, target_path)
        print(f"   renamed 已重命名: {entry} → {target_name}")
        renamed += 1

    if candidates == 0:
        print("   所有树文件已是 .nwk / All tree files already use .nwk")
    else:
        print(f"   共重命名 {renamed} 个文件 / Renamed {renamed} file(s)")
    print()

    return renamed


def collect_files(directory):
    fasta_map = {}
    tree_map = {}
    fasta_dups = defaultdict(list)
    tree_dups = defaultdict(list)
    note_map = defaultdict(list)
    warning_messages = []

    for entry in sorted(os.listdir(directory)):
        full_path = os.path.join(directory, entry)
        if not os.path.isfile(full_path):
            continue

        stem, ext = os.path.splitext(entry)
        ext = ext.lower()
        if ext in FASTA_EXTS:
            fasta_dups[stem].append(full_path)
            fasta_map.setdefault(stem, full_path)
        elif ext in TREE_EXTS:
            tree_dups[stem].append(full_path)
            tree_map.setdefault(stem, full_path)

    for stem in sorted(fasta_dups):
        paths = fasta_dups[stem]
        if len(paths) > 1:
            basenames = [os.path.basename(path) for path in paths]
            note_map[stem].append(f"multiple_fasta={'|'.join(basenames)}")
            used = os.path.basename(fasta_map[stem])
            ignored = ", ".join(name for name in basenames if name != used)
            warning_messages.append(
                f"⚠️  {stem:<25s} MULTIPLE_FASTA: using 使用 {used}, ignoring 忽略 {ignored}"
            )

    for stem in sorted(tree_dups):
        paths = tree_dups[stem]
        if len(paths) > 1:
            basenames = [os.path.basename(path) for path in paths]
            note_map[stem].append(f"multiple_tree={'|'.join(basenames)}")
            used = os.path.basename(tree_map[stem])
            ignored = ", ".join(name for name in basenames if name != used)
            warning_messages.append(
                f"⚠️  {stem:<25s} MULTIPLE_TREE:  using 使用 {used}, ignoring 忽略 {ignored}"
            )

    return FileCollection(
        fasta_map=fasta_map,
        tree_map=tree_map,
        note_map=note_map,
        warning_messages=warning_messages,
    )


def _print_collection_warnings(collection):
    for message in collection.warning_messages:
        print(message)
    if collection.warning_messages:
        print()


def _write_mismatch_log(log_path, prefix, fasta_path, tree_path, fa_taxa, tr_taxa, tree_only, fasta_only):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"PREFIX: {prefix}\n")
        handle.write(f"FASTA:  {fasta_path}\n")
        handle.write(f"TREE:   {tree_path}\n")
        handle.write(f"Counts: TREE={len(tr_taxa)} FASTA={len(fa_taxa)}\n\n")
        handle.write(
            f"In TREE but NOT in FASTA / 在树中但不在FASTA中 ({len(tree_only)}):\n"
        )
        for taxon in tree_only:
            handle.write(f"  {taxon}\n")
        handle.write(
            f"\nIn FASTA but NOT in TREE / 在FASTA中但不在树中 ({len(fasta_only)}):\n"
        )
        for taxon in fasta_only:
            handle.write(f"  {taxon}\n")
        handle.write(
            "\nHint 提示: also check for case / underscore / dot / dash differences "
            "if counts look close.\n"
        )
        handle.write(
            "如果数量接近，请检查大小写、下划线、点、横线等格式差异。\n"
        )
        handle.write(
            "\nNote 注意: FASTA taxa = first token after '>' (space-delimited); "
            "tree may use quoted names with spaces.\n"
        )
        handle.write(
            "FASTA 物种名 = '>' 后第一个空白前的部分；树中引号包裹的名称会保留空格。\n"
        )


def run_check_directory(directory, logdir, do_clean=False, do_normalize_ext=False, prefix_dir=""):
    if do_normalize_ext:
        normalize_tree_extensions(directory)

    collection = collect_files(directory)
    _print_collection_warnings(collection)

    if do_clean:
        print("🔧 清理树文件中的内部节点标签... / Cleaning internal node labels...")
        cleaned = 0
        for prefix in sorted(collection.tree_map):
            tree_path = collection.tree_map[prefix]
            backup, removed = clean_internal_labels(tree_path)
            if removed > 0:
                print(f"   cleaned 已清理: {tree_path} (backup 备份: {backup})")
                cleaned += 1
        if cleaned == 0:
            print("   未发现需要清理的内部节点标签 / No internal node labels found")
        else:
            print(f"   共清理 {cleaned} 个树文件 / Cleaned {cleaned} tree file(s)")
        print()

    rows = []
    matched = 0
    ok = 0
    bad = 0
    prefixes = sorted(set(collection.fasta_map) | set(collection.tree_map))

    if not prefixes:
        print(f"⚠️  目录为空或无可识别文件 / No FASTA or tree files found: {directory}")
        print("----")
        print("Paired 已配对: 0   OK 通过: 0   MISMATCH 不一致: 0")
        return rows, 0

    for prefix in prefixes:
        display_prefix = _format_prefix(prefix_dir, prefix)
        notes = ";".join(collection.note_map.get(prefix, []))
        fasta_path = collection.fasta_map.get(prefix, "")
        tree_path = collection.tree_map.get(prefix, "")

        if fasta_path and tree_path:
            matched += 1
            fa_taxa = set(extract_fasta_taxa(fasta_path))
            tr_taxa = set(extract_tree_taxa(tree_path))
            tree_only = sorted(tr_taxa - fa_taxa)
            fasta_only = sorted(fa_taxa - tr_taxa)

            if not tree_only and not fasta_only:
                print(f"✅ {prefix:<25s} OK (n={len(tr_taxa)})")
                ok += 1
                status = "OK"
                log_path = ""
            else:
                bad += 1
                log_path = os.path.join(logdir, f"{prefix}.mismatch.txt")
                _write_mismatch_log(
                    log_path,
                    prefix,
                    fasta_path,
                    tree_path,
                    fa_taxa,
                    tr_taxa,
                    tree_only,
                    fasta_only,
                )
                print(
                    f"❌ {prefix:<25s} MISMATCH 不一致 "
                    f"(tree_only={len(tree_only)} fasta_only={len(fasta_only)}) -> {log_path}"
                )
                status = "MISMATCH"

            rows.append(
                {
                    "prefix": display_prefix,
                    "status": status,
                    "notes": notes,
                    "n_fasta": len(fa_taxa),
                    "n_tree": len(tr_taxa),
                    "n_tree_only": len(tree_only),
                    "n_fasta_only": len(fasta_only),
                    "fasta_path": fasta_path,
                    "tree_path": tree_path,
                    "log_path": log_path,
                }
            )
            continue

        if fasta_path:
            fa_taxa = extract_fasta_taxa(fasta_path)
            print(f"⚠️  {prefix:<25s} NO_TREE 缺树文件  (have FASTA 有FASTA: {fasta_path})")
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "NO_TREE",
                    "notes": notes,
                    "n_fasta": len(fa_taxa),
                    "n_tree": 0,
                    "n_tree_only": 0,
                    "n_fasta_only": len(fa_taxa),
                    "fasta_path": fasta_path,
                    "tree_path": "",
                    "log_path": "",
                }
            )
        else:
            tr_taxa = extract_tree_taxa(tree_path)
            print(f"⚠️  {prefix:<25s} NO_FASTA 缺FASTA (have TREE 有树文件: {tree_path})")
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "NO_FASTA",
                    "notes": notes,
                    "n_fasta": 0,
                    "n_tree": len(tr_taxa),
                    "n_tree_only": len(tr_taxa),
                    "n_fasta_only": 0,
                    "fasta_path": "",
                    "tree_path": tree_path,
                    "log_path": "",
                }
            )

    print("----")
    print(f"Paired 已配对: {matched}   OK 通过: {ok}   MISMATCH 不一致: {bad}")
    if bad > 0:
        print(f"Mismatch logs 不一致日志: {logdir}/")

    return rows, bad


def run_clean_directory(directory, prefix_dir=""):
    collection = collect_files(directory)
    _print_collection_warnings(collection)

    print("🔧 清理树文件中的内部节点标签... / Cleaning internal node labels...")
    rows = []
    cleaned = 0
    no_tree = 0
    prefixes = sorted(set(collection.fasta_map) | set(collection.tree_map))

    if not prefixes:
        print("   未发现树文件 / No tree files found")
        return rows, 0

    for prefix in prefixes:
        display_prefix = _format_prefix(prefix_dir, prefix)
        tree_path = collection.tree_map.get(prefix, "")
        if not tree_path:
            fasta_path = collection.fasta_map.get(prefix, "")
            print(f"⚠️  {prefix:<25s} NO_TREE 缺树文件  (have FASTA 有FASTA: {fasta_path})")
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "NO_TREE",
                    "n_labels_removed": 0,
                    "tree_path": "",
                    "backup_path": "",
                }
            )
            no_tree += 1
            continue

        backup, removed = clean_internal_labels(tree_path)
        if removed > 0:
            print(f"   cleaned 已清理: {tree_path} (backup 备份: {backup})")
            cleaned += 1
            status = "CLEANED"
        else:
            status = "NO_LABELS"
            backup = ""

        rows.append(
            {
                "prefix": display_prefix,
                "status": status,
                "n_labels_removed": removed,
                "tree_path": tree_path,
                "backup_path": backup,
            }
        )

    if cleaned == 0:
        print("   未发现需要清理的内部节点标签 / No internal node labels found")
    else:
        print(f"   共清理 {cleaned} 个树文件 / Cleaned {cleaned} tree file(s)")

    return rows, no_tree


def _quote_tree_label(label, preferred_quote=None):
    if preferred_quote:
        quote = preferred_quote
    elif "'" not in label:
        quote = "'"
    else:
        quote = '"'

    if quote == "'":
        escaped = label.replace("'", "''")
    else:
        escaped = label.replace('"', '""')
    return f"{quote}{escaped}{quote}"


def _render_tree_label(label, quote_char=None):
    if quote_char:
        return _quote_tree_label(label, quote_char)
    if UNQUOTED_TREE_LABEL_RE.match(label):
        return label
    return _quote_tree_label(label)


def rewrite_tree_leaves(text, rename_map):
    out = []
    i = 0
    n = len(text)
    expect_leaf = True

    while i < n:
        ch = text[i]

        if ch == "(":
            expect_leaf = True
            out.append(ch)
            i += 1
        elif ch == ")":
            expect_leaf = False
            out.append(ch)
            i += 1
        elif ch == ",":
            expect_leaf = True
            out.append(ch)
            i += 1
        elif ch == ";":
            out.append(ch)
            i += 1
        elif ch == "[":
            start = i
            depth = 1
            i += 1
            while i < n and depth > 0:
                if text[i] == "[":
                    depth += 1
                elif text[i] == "]":
                    depth -= 1
                i += 1
            out.append(text[start:i])
        elif ch == ":":
            start = i
            i += 1
            while i < n and text[i] in "0123456789.eE+-":
                i += 1
            out.append(text[start:i])
        elif ch == "'":
            label, i = _consume_quoted_label(text, i, "'")
            if expect_leaf:
                label = rename_map.get(label, label)
                out.append(_render_tree_label(label, "'"))
            else:
                out.append(_render_tree_label(label, "'"))
        elif ch == '"':
            label, i = _consume_quoted_label(text, i, '"')
            if expect_leaf:
                label = rename_map.get(label, label)
                out.append(_render_tree_label(label, '"'))
            else:
                out.append(_render_tree_label(label, '"'))
        elif ch in " \t\r\n":
            out.append(ch)
            i += 1
        else:
            start = i
            while i < n and text[i] not in ":,();[] \t\r\n":
                i += 1
            label = text[start:i]
            if expect_leaf:
                out.append(_render_tree_label(rename_map.get(label, label)))
            else:
                out.append(label)

    return "".join(out)


def rewrite_fasta_headers(text, rename_map):
    out_lines = []
    for line in text.splitlines():
        if not line.startswith(">"):
            out_lines.append(line)
            continue

        body = line[1:]
        stripped = body.lstrip()
        if not stripped:
            out_lines.append(line)
            continue

        leading_spaces = len(body) - len(stripped)
        token, _, rest = stripped.partition(" ")
        new_token = rename_map.get(token, token)
        if rest:
            new_line = ">" + (" " * leading_spaces) + new_token + " " + rest
        else:
            new_line = ">" + (" " * leading_spaces) + new_token
        out_lines.append(new_line)

    trailing_newline = "\n" if text.endswith("\n") else ""
    return "\n".join(out_lines) + trailing_newline


def load_mapping_file(filepath):
    with open(filepath, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["old_name", "new_name"]:
            raise ValueError("mapping TSV must contain header: old_name<TAB>new_name")
        mapping = {}
        for row in reader:
            old_name = (row.get("old_name") or "").strip()
            new_name = (row.get("new_name") or "").strip()
            if not old_name:
                continue
            mapping[old_name] = new_name
    return mapping


def parse_abbrev_parts(spec):
    parts = []
    for raw in spec.split(","):
        token = raw.strip()
        if not token:
            raise ValueError("empty token in --parts")
        if token.lower() == "full":
            parts.append("full")
            continue
        value = int(token)
        if value < 0:
            raise ValueError("--parts must use non-negative integers or 'full'")
        parts.append(value)
    return parts


def transform_name(name, args):
    if args.rule == "upper":
        return name.upper()
    if args.rule == "lower":
        return name.lower()
    if args.rule == "title":
        return name[:1].upper() + name[1:].lower()
    if args.rule == "abbrev":
        tokens = [token for token in re.split(r"[_. ]+", name.strip()) if token]
        pieces = []
        for token, part_spec in zip(tokens, args.parsed_parts):
            if part_spec == "full":
                pieces.append(token)
            else:
                pieces.append(token[:part_spec])
        return args.sep.join(pieces)
    if args.rule == "truncate-dot":
        return ".".join(name.split(".")[:args.keep_parts])
    if args.rule == "join-dot":
        return name.replace(".", "")
    if args.rule == "replace":
        return name.replace(args.replace_from, args.replace_to)
    if args.rule == "strip-version":
        return re.sub(r"\.\d+$", "", name)
    if args.rule == "from-file":
        if name in args.name_mapping:
            return args.name_mapping[name]
        if args.strict_mapping:
            raise KeyError(name)
        return name
    raise ValueError(f"unsupported rule: {args.rule}")


def _write_mapping_tsv(mapping_path, rename_map):
    rows = [{"old_name": old_name, "new_name": rename_map[old_name]} for old_name in sorted(rename_map)]
    _write_tsv(mapping_path, ["old_name", "new_name"], rows)


def _detect_rename_collisions(rename_map):
    collisions = defaultdict(list)
    for old_name, new_name in rename_map.items():
        collisions[new_name].append(old_name)
    return {
        new_name: sorted(old_names)
        for new_name, old_names in collisions.items()
        if len(old_names) > 1
    }


def run_rename_directory(directory, args, prefix_dir=""):
    collection = collect_files(directory)
    _print_collection_warnings(collection)

    rows = []
    failure_units = 0
    prefixes = sorted(set(collection.fasta_map) | set(collection.tree_map))

    if not prefixes:
        print(f"⚠️  目录为空或无可识别文件 / No FASTA or tree files found: {directory}")
        return rows, 0

    for prefix in prefixes:
        display_prefix = _format_prefix(prefix_dir, prefix)
        fasta_path = collection.fasta_map.get(prefix, "")
        tree_path = collection.tree_map.get(prefix, "")
        mapping_path = os.path.join(directory, f"{prefix}_rename_mapping.tsv")

        if not tree_path:
            print(f"⚠️  {prefix:<25s} NO_TREE 缺树文件  (have FASTA 有FASTA: {fasta_path})")
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "NO_TREE",
                    "n_renamed": 0,
                    "n_unchanged": 0,
                    "n_collisions": 0,
                    "mapping_path": "",
                }
            )
            continue

        if not fasta_path:
            print(f"⚠️  {prefix:<25s} NO_FASTA 缺FASTA (have TREE 有树文件: {tree_path})")
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "NO_FASTA",
                    "n_renamed": 0,
                    "n_unchanged": 0,
                    "n_collisions": 0,
                    "mapping_path": "",
                }
            )
            continue

        fasta_taxa = set(extract_fasta_taxa(fasta_path))
        tree_taxa = set(extract_tree_taxa(tree_path))
        all_taxa = sorted(fasta_taxa | tree_taxa)

        if args.rule == "from-file":
            unused = sorted(set(args.name_mapping) - set(all_taxa))
            if unused:
                print(
                    f"⚠️  {prefix:<25s} UNUSED_MAPPING 未使用映射: "
                    f"{', '.join(unused)}"
                )

        rename_map = {}
        missing = []
        for old_name in all_taxa:
            try:
                rename_map[old_name] = transform_name(old_name, args)
            except KeyError:
                missing.append(old_name)

        if missing:
            print(
                f"❌ {prefix:<25s} STRICT_MAPPING 缺少映射: "
                f"{', '.join(sorted(missing))}"
            )
            print("Aborted. No files were modified.")
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "ERROR",
                    "n_renamed": 0,
                    "n_unchanged": len(all_taxa),
                    "n_collisions": 0,
                    "mapping_path": "",
                }
            )
            failure_units += 2
            continue

        collisions = _detect_rename_collisions(rename_map)

        if collisions:
            print("❌ Rename collision detected / 检测到重命名冲突:")
            for new_name in sorted(collisions):
                print(f"   {new_name}  ← {', '.join(collisions[new_name])}")
            print("Aborted. No files were modified.")
            print(
                "建议: use --rule from-file with manual mapping to resolve / "
                "使用 --rule from-file 手动指定 mapping"
            )
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "COLLISION",
                    "n_renamed": 0,
                    "n_unchanged": len(all_taxa),
                    "n_collisions": len(collisions),
                    "mapping_path": "",
                }
            )
            failure_units += 2
            continue

        n_renamed = sum(1 for old_name, new_name in rename_map.items() if old_name != new_name)
        n_unchanged = len(rename_map) - n_renamed
        _write_mapping_tsv(mapping_path, rename_map)

        if n_renamed == 0:
            print(f"✅ {prefix:<25s} NO_CHANGE 无实际改动 -> {mapping_path}")
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "NO_CHANGE",
                    "n_renamed": 0,
                    "n_unchanged": n_unchanged,
                    "n_collisions": 0,
                    "mapping_path": mapping_path,
                }
            )
            continue

        if args.dry_run:
            print(
                f"👀 {prefix:<25s} DRY_RUN 预览改名 "
                f"(renamed={n_renamed} unchanged={n_unchanged}) -> {mapping_path}"
            )
            rows.append(
                {
                    "prefix": display_prefix,
                    "status": "DRY_RUN",
                    "n_renamed": n_renamed,
                    "n_unchanged": n_unchanged,
                    "n_collisions": 0,
                    "mapping_path": mapping_path,
                }
            )
            continue

        fasta_text = _read_text(fasta_path)
        tree_text = _read_text(tree_path)
        new_fasta_text = rewrite_fasta_headers(fasta_text, rename_map)
        new_tree_text = rewrite_tree_leaves(tree_text, rename_map)

        _safe_backup(fasta_path)
        _safe_backup(tree_path)
        _write_text(fasta_path, new_fasta_text)
        _write_text(tree_path, new_tree_text)

        print(
            f"✏️  {prefix:<25s} APPLIED 已写入 "
            f"(renamed={n_renamed} unchanged={n_unchanged}) -> {mapping_path}"
        )
        rows.append(
            {
                "prefix": display_prefix,
                "status": "APPLIED",
                "n_renamed": n_renamed,
                "n_unchanged": n_unchanged,
                "n_collisions": 0,
                "mapping_path": mapping_path,
            }
        )

    return rows, failure_units


def iter_target_directories(directory, batch, logdir_name):
    if not batch:
        return [("", directory)]

    subdirs = sorted(
        entry
        for entry in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, entry)) and entry != logdir_name
    )
    return [(subdir, os.path.join(directory, subdir)) for subdir in subdirs]


def resolve_logdir(directory, logdir):
    if os.path.isabs(logdir):
        return logdir
    return os.path.join(directory, logdir)


def run_check(args):
    targets = iter_target_directories(args.directory, args.batch, args.logdir)
    if args.batch and not targets:
        print(f"未找到子目录 / No subdirectories found: {args.directory}")
        return 1

    all_rows = []
    total_failures = 0
    for subdir, directory in targets:
        prefix_dir = subdir if args.batch else ""
        if args.batch:
            print(f"===== {subdir} =====")
        logdir = resolve_logdir(directory, args.logdir)
        rows, failures = run_check_directory(
            directory,
            logdir,
            do_clean=args.clean,
            do_normalize_ext=args.normalize_ext,
            prefix_dir=prefix_dir,
        )
        all_rows.extend(rows)
        total_failures += failures
        if args.batch:
            print()

    if args.report:
        _write_tsv(args.report, CHECK_REPORT_FIELDS, all_rows)

    return min(total_failures, 255)


def run_clean(args):
    targets = iter_target_directories(args.directory, args.batch, args.logdir)
    if args.batch and not targets:
        print(f"未找到子目录 / No subdirectories found: {args.directory}")
        return 1

    all_rows = []
    total_failures = 0
    for subdir, directory in targets:
        prefix_dir = subdir if args.batch else ""
        if args.batch:
            print(f"===== {subdir} =====")
        rows, failures = run_clean_directory(directory, prefix_dir=prefix_dir)
        all_rows.extend(rows)
        total_failures += failures
        if args.batch:
            print()

    if args.report:
        _write_tsv(args.report, CLEAN_REPORT_FIELDS, all_rows)

    return min(total_failures, 255)


def run_rename(args):
    targets = iter_target_directories(args.directory, args.batch, args.logdir)
    if args.batch and not targets:
        print(f"未找到子目录 / No subdirectories found: {args.directory}")
        return 1

    all_rows = []
    total_failures = 0
    for subdir, directory in targets:
        prefix_dir = subdir if args.batch else ""
        if args.batch:
            print(f"===== {subdir} =====")
        rows, failures = run_rename_directory(directory, args, prefix_dir=prefix_dir)
        all_rows.extend(rows)
        total_failures += failures
        if args.batch:
            print()

    if args.report:
        _write_tsv(args.report, RENAME_REPORT_FIELDS, all_rows)

    return min(total_failures, 255)


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "批量检查 FASTA 与树文件之间的物种名一致性，并提供 clean / rename 子命令\n"
            "Batch check taxon name consistency between FASTA and tree files"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="待处理目录 / Directory to process (默认/default: current dir)",
    )
    common.add_argument(
        "--batch",
        action="store_true",
        help="批量模式: 遍历第一层子目录 / Batch mode: iterate first-level subdirs",
    )
    common.add_argument(
        "--logdir",
        default="taxa_mismatch_logs",
        help="不一致日志目录 / Mismatch log directory (默认/default: taxa_mismatch_logs)",
    )
    common.add_argument(
        "--report",
        help="输出 TSV 报告 / Write TSV report",
    )

    check_parser = subparsers.add_parser(
        "check",
        parents=[common],
        help="检查 FASTA 与树叶节点的一致性 / Check FASTA vs tree taxa",
    )
    check_parser.add_argument(
        "--clean",
        action="store_true",
        help="检查前先清理树内部节点标签 / Clean internal labels before checking",
    )
    check_parser.add_argument(
        "--normalize-ext",
        action="store_true",
        help="检查前将 .newick/.treefile/.tre 树文件扩展名统一为 .nwk / Normalize tree extensions to .nwk before checking",
    )

    subparsers.add_parser(
        "clean",
        parents=[common],
        help="仅清理树内部节点标签 / Clean internal node labels only",
    )

    rename_parser = subparsers.add_parser(
        "rename",
        parents=[common],
        help="同步重命名 FASTA headers + 树叶节点 / Rename FASTA headers and tree leaves",
    )
    rename_parser.add_argument(
        "--rule",
        required=True,
        help=(
            "重命名规则 / Rename rule: upper, lower, title(first-char only), "
            "abbrev, truncate-dot, join-dot, replace, strip-version, from-file"
        ),
    )
    rename_parser.add_argument("--dry-run", action="store_true", help="只预览，不改文件 / Preview only")
    rename_parser.add_argument("--parts", help="abbrev 规则: 截断长度列表 / abbrev parts")
    rename_parser.add_argument(
        "--sep",
        default="_",
        help='abbrev 规则: 输出分隔符 / abbrev separator (use "" for no separator)',
    )
    rename_parser.add_argument(
        "--keep-parts",
        type=int,
        help="truncate-dot 规则: 保留前 N 段 / truncate-dot keep N parts",
    )
    rename_parser.add_argument("--from", dest="replace_from", help="replace 规则: old string")
    rename_parser.add_argument("--to", dest="replace_to", default="", help="replace 规则: new string")
    rename_parser.add_argument("--map", help="from-file 规则: TSV mapping 文件 / TSV mapping file")
    rename_parser.add_argument(
        "--strict-mapping",
        action="store_true",
        help="from-file 规则: 缺失映射时报错 / fail when mapping is incomplete",
    )

    return parser


def normalize_argv(argv):
    if not argv:
        return ["check"]
    if argv[0] in {"check", "clean", "rename", "-h", "--help"}:
        return argv
    return ["check", *argv]


def validate_args(args):
    if args.command != "rename":
        return

    valid_rules = {
        "upper",
        "lower",
        "title",
        "abbrev",
        "truncate-dot",
        "join-dot",
        "replace",
        "strip-version",
        "from-file",
    }
    if args.rule not in valid_rules:
        raise ValueError(f"unsupported --rule: {args.rule}")

    if args.rule == "abbrev":
        if not args.parts:
            raise ValueError("--rule abbrev requires --parts")
        args.parsed_parts = parse_abbrev_parts(args.parts)
    else:
        args.parsed_parts = []

    if args.rule == "truncate-dot":
        if args.keep_parts is None or args.keep_parts < 0:
            raise ValueError("--rule truncate-dot requires --keep-parts N (N >= 0)")

    if args.rule == "replace":
        if args.replace_from is None:
            raise ValueError("--rule replace requires --from STR")
        if args.replace_from == "":
            raise ValueError("--from cannot be empty")

    if args.rule == "from-file":
        if not args.map:
            raise ValueError("--rule from-file requires --map FILE")
        args.name_mapping = load_mapping_file(args.map)
    else:
        args.name_mapping = {}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = build_parser()

    try:
        args = parser.parse_args(normalize_argv(argv))
        if args.command is None:
            parser.print_help()
            return 0
        validate_args(args)
    except Exception as exc:
        print(f"Error / 错误: {exc}", file=sys.stderr)
        return 2

    if args.command == "check":
        return run_check(args)
    if args.command == "clean":
        return run_clean(args)
    if args.command == "rename":
        return run_rename(args)
    parser.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
