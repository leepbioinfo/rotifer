import re
from pathlib import Path

import pandas as pd
from rotifer.devel.beta.sequence import sequence


def _parse_int_list(x):
    return [int(n) for n in re.findall(r"\d+", str(x))]


def _decompress_dssp(x):
    if pd.isna(x):
        return ""

    parts = re.findall(r"(\d+)([A-Za-z-])", str(x))
    return "".join(char * int(n) for n, char in parts)


def _format_unaligned(text, mode="lower"):
    if mode == "keep":
        return text
    if mode == "lower":
        return text.lower()
    if mode == "gap":
        return "-" * len(text)
    if mode == "remove":
        return ""

    raise ValueError("unaligned must be one of: 'keep', 'lower', 'gap', 'remove'")


def _read_dali_tsv(path, max_hits=None):
    path = Path(path)

    query_seq = None
    query_dssp = None

    with path.open() as f:
        for line in f:
            if line.startswith("# SEQUENCE:"):
                query_seq = line.split(":", 1)[1].strip()
            elif line.startswith("# DSSP:"):
                query_dssp = line.split(":", 1)[1].strip()

    if query_seq is None:
        raise ValueError("Could not find '# SEQUENCE:' line in DALI file.")

    dali = pd.read_csv(path, sep="\t", comment="#")

    if max_hits is not None:
        dali = dali.head(max_hits).copy()

    return dali, query_seq, query_dssp


def _row_blocks(row):
    """
    DALI TSV qstarts/sstarts are 0-based.
    Convert to 1-based coordinates for the slicing logic below.
    """
    qstarts = [x + 1 for x in _parse_int_list(row["qstarts"])]
    sstarts = [x + 1 for x in _parse_int_list(row["sstarts"])]
    lengths = _parse_int_list(row["lengths"])

    return list(zip(qstarts, sstarts, lengths))


def _collect_insert_sizes(dali, query_len):
    insert_sizes = [0] * (query_len + 1)

    for _, row in dali.iterrows():
        subject_seq = row["sbjct-sequence"]
        prev_s_end = 1

        for q_start, s_start, length in _row_blocks(row):
            q_insert_pos = q_start - 1
            subject_insertion_len = s_start - prev_s_end

            insert_sizes[q_insert_pos] = max(
                insert_sizes[q_insert_pos],
                subject_insertion_len,
            )

            prev_s_end = s_start + length

        tail_len = len(subject_seq) - prev_s_end + 1
        insert_sizes[query_len] = max(
            insert_sizes[query_len],
            max(0, tail_len),
        )

    return insert_sizes


def _make_star_query(query_seq, insert_sizes):
    aln = []

    for i, char in enumerate(query_seq):
        aln.append("-" * insert_sizes[i])
        aln.append(char.upper())

    aln.append("-" * insert_sizes[len(query_seq)])

    return "".join(aln)


def _make_star_subject(row, query_seq, insert_sizes, unaligned="lower"):
    subject_seq = row["sbjct-sequence"]
    query_len = len(query_seq)

    subject_on_query = ["-"] * query_len
    subject_inserts = [""] * (query_len + 1)

    prev_s_end = 1

    for q_start, s_start, length in _row_blocks(row):
        q_insert_pos = q_start - 1

        insertion = subject_seq[prev_s_end - 1:s_start - 1]
        subject_inserts[q_insert_pos] += _format_unaligned(
            insertion,
            mode=unaligned,
        )

        for offset in range(length):
            q_idx = q_start - 1 + offset
            s_idx = s_start - 1 + offset

            if 0 <= q_idx < query_len and 0 <= s_idx < len(subject_seq):
                subject_on_query[q_idx] = subject_seq[s_idx].upper()

        prev_s_end = s_start + length

    tail = subject_seq[prev_s_end - 1:]
    subject_inserts[query_len] += _format_unaligned(tail, mode=unaligned)

    aln = []

    for i, char in enumerate(subject_on_query):
        insertion = subject_inserts[i]
        padded_insertion = insertion + "-" * (insert_sizes[i] - len(insertion))

        aln.append(padded_insertion)
        aln.append(char)

    tail = subject_inserts[query_len]
    padded_tail = tail + "-" * (insert_sizes[query_len] - len(tail))
    aln.append(padded_tail)

    return "".join(aln)


def _make_star_annotation(
    row,
    query_seq,
    insert_sizes,
    annotation_seq,
    unaligned="lower",
):
    query_len = len(query_seq)

    ann_on_query = ["-"] * query_len
    ann_inserts = [""] * (query_len + 1)

    prev_s_end = 1

    for q_start, s_start, length in _row_blocks(row):
        q_insert_pos = q_start - 1

        insertion = annotation_seq[prev_s_end - 1:s_start - 1]
        ann_inserts[q_insert_pos] += _format_unaligned(
            insertion,
            mode=unaligned,
        )

        for offset in range(length):
            q_idx = q_start - 1 + offset
            s_idx = s_start - 1 + offset

            if 0 <= q_idx < query_len and 0 <= s_idx < len(annotation_seq):
                ann_on_query[q_idx] = annotation_seq[s_idx].upper()

        prev_s_end = s_start + length

    tail = annotation_seq[prev_s_end - 1:]
    ann_inserts[query_len] += _format_unaligned(tail, mode=unaligned)

    aln = []

    for i, char in enumerate(ann_on_query):
        insertion = ann_inserts[i]
        padded_insertion = insertion + "-" * (insert_sizes[i] - len(insertion))

        aln.append(padded_insertion)
        aln.append(char)

    tail = ann_inserts[query_len]
    padded_tail = tail + "-" * (insert_sizes[query_len] - len(tail))
    aln.append(padded_tail)

    return "".join(aln)


def _get_subject_description(row):
    if "description" in row.index and pd.notna(row["description"]):
        return row["description"]

    if "sbjct-description" in row.index and pd.notna(row["sbjct-description"]):
        return row["sbjct-description"]

    return row["sbjct"]


def dali_tsv(
    path,
    max_hits=None,
    n_structures=None,
    add_dssp=True,
    unaligned="lower",
):
    dali, query_seq, query_dssp = _read_dali_tsv(path, max_hits=max_hits)

    if dali.empty:
        raise ValueError("No DALI hits available after filtering.")

    # sequences included in the alignment
    dali_seq = dali.copy()

    # structures with DSSP annotation rows
    if n_structures is not None:
        dali_dssp = dali.head(n_structures).copy()
    else:
        dali_dssp = dali.copy()

    if unaligned == "remove":
        insert_sizes = [0] * (len(query_seq) + 1)
    else:
        insert_sizes = _collect_insert_sizes(dali_seq, len(query_seq))

    query_id = dali_seq["query"].iloc[0]

    query_description = (
        dali_seq["query-description"].iloc[0]
        if "query-description" in dali_seq.columns
        else query_id
    )

    query_aln = _make_star_query(query_seq, insert_sizes)

    annotation_rows = []
    sequence_rows = []

    if add_dssp and query_dssp is not None:
        annotation_rows.append({
            "id": f"ss_from:{query_id}",
            "sequence": _make_star_query(query_dssp, insert_sizes),
            "type": "residue_annotation",
            "description": query_description,
        })

    sequence_rows.append({
        "id": query_id,
        "sequence": query_aln,
        "type": "sequence",
        "description": query_description,
    })

    # output all sequence rows from max_hits
    for _, row in dali_seq.iterrows():
        subject_id = row["sbjct"]
        subject_description = _get_subject_description(row)

        subject_aln = _make_star_subject(
            row=row,
            query_seq=query_seq,
            insert_sizes=insert_sizes,
            unaligned=unaligned,
        )

        sequence_rows.append({
            "id": subject_id,
            "sequence": subject_aln,
            "type": "sequence",
            "description": subject_description,
        })

    # output DSSP only for n_structures
    if add_dssp:
        for _, row in dali_dssp.iterrows():
            subject_id = row["sbjct"]
            subject_description = _get_subject_description(row)

            subject_dssp = _decompress_dssp(row["compressed-sbjct-dssp"])

            annotation_rows.append({
                "id": f"ss_from:{subject_id}",
                "sequence": _make_star_annotation(
                    row=row,
                    query_seq=query_seq,
                    insert_sizes=insert_sizes,
                    annotation_seq=subject_dssp,
                    unaligned=unaligned,
                ),
                "type": "residue_annotation",
                "description": subject_description,
            })

    df = pd.DataFrame(annotation_rows + sequence_rows)

    aln = sequence()
    aln.df = df
    aln._reset()

    return aln