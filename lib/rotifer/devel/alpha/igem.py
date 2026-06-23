"""
operon_fig.py
=============

Draw gene-neighborhood ("operon") figures from a long-format table of
genes/proteins that have already been grouped into genomic blocks (e.g.
one block per hit of interest plus its surrounding genomic context).

Each block is drawn as one row of the figure:

    [ text label ]  [gene] [gene] [QUERY] [gene] [gene] ...

The text label is a 3-line box with the query protein id, the block id
and the organism name. Genes are drawn as right-/left-pointing arrows
depending on strand and colored by a domain/annotation label (Pfam,
Aravind, or free-text product, depending on `label_col`).

Required input columns
-----------------------
pid     : protein id (used in the row label and to find the query gene)
strand  : +1 / -1 (controls arrow direction and orientation normalization)
query   : 1/'1'/True marks the query gene of a block; everything else is
          treated as a non-query gene. Optional -- if absent, no gene is
          treated as a query.
... plus whatever `group_col`, `org_col` and `label_col` point at.

Public entry point
-------------------
operon_fig_patched(df, ...) -> pandas.DataFrame
    Builds the figure, writes it to `output_file`, and returns the
    working copy of `df` that was used to build it (useful for
    debugging/inspection). Kept under its original name so existing
    calling code does not need to change.

Everything else below is a small, independently testable helper
function. None of them are defined *inside* another function (the
original had a `pad_and_escape` closure) -- they all live at module
level so they can be imported and unit tested on their own, and so a
future change to one of them does not require re-reading the whole
pipeline.

    resolve_domain_labels        fill the 'domain' column (labels/colors)
    flag_query_rows              add the boolean 'is_query' column
    prepare_dataframe            run the two helpers above + housekeeping
    compute_label_width          shared padding width for row labels
    pad_and_escape               pad + HTML-escape one label string
    build_row_label_html         the 3-line HTML label for one block
    build_color_map              domain -> fill color, incl. user overrides
    normalize_block_strand       optionally mirror a block to a common
                                  query orientation
    collapsed_triangle_orientation  Graphviz 'orientation' degrees for a
                                  collapsed opposite-strand gene
    gene_node_style               Graphviz node attributes for one gene
    add_block_to_graph            add one full row (label + genes) to the
                                  graph, with optional left/right padding
    chain_align_nodes             pull one node per row into the same
                                  visual column via high-weight invisible
                                  edges
    operon_fig_patched            main orchestrator (see above)
"""

import html

import numpy as np
import pandas as pd
import pygraphviz as pgv
import seaborn as sns

# Domain/annotation values that never get an automatic color, because
# they are generic/structural rather than informative (signal peptides,
# transmembrane regions, etc.) or simply mean "no annotation".
DEFAULT_IGNORE_DOMAINS = ['TM', 'SP', 'LP', 'LIPO', 'SIG']


# ---------------------------------------------------------------------------
# Column preparation
# ---------------------------------------------------------------------------

def resolve_domain_labels(df, label_col='pfam'):
    """
    Build the 'domain' column that gene node labels/colors are based on.

    When `label_col == 'pfam'` (the default), a gene with no Pfam hit
    falls back to its 'aravind' annotation, then to its free-text
    'product' description, and finally to the literal string 'unk' if
    none of those are available either. This mirrors the common
    situation where only some genes in a neighborhood have a Pfam domain.

    For any other `label_col`, that column is used as-is (missing values
    become 'unk'); if the column does not exist at all, every gene gets
    'unk'.

    Parameters
    ----------
    df : pandas.DataFrame
    label_col : str

    Returns
    -------
    pandas.DataFrame
        Copy of `df` with a 'domain' column added.
    """
    out = df.copy()

    if label_col == 'pfam':
        domain = out['pfam'] if 'pfam' in out.columns else pd.Series(np.nan, index=out.index)
        if 'aravind' in out.columns:
            domain = domain.fillna(out['aravind'])
        if 'product' in out.columns:
            domain = domain.fillna(out['product'])
        out['domain'] = domain.fillna('unk')
    elif label_col in out.columns:
        out['domain'] = out[label_col].fillna('unk')
    else:
        out['domain'] = 'unk'

    return out


def flag_query_rows(df):
    """
    Add a boolean 'is_query' column derived from the raw 'query' column.

    1, '1' and True are all treated as "this is the query gene"; a
    missing 'query' column means no gene in `df` is a query.
    """
    out = df.copy()
    if 'query' in out.columns:
        out['is_query'] = out['query'].isin([1, '1', True])
    else:
        out['is_query'] = False
    return out


def prepare_dataframe(df, group_col='block_id', org_col='organism', label_col='pfam'):
    """
    Normalize a raw input table into the columns the rest of this module
    relies on: 'ID' (block id), 'org_name', 'domain', 'is_query' and
    'pid_order' (an integer 0..n_blocks-1, in first-seen order, used to
    group rows into figure rows).

    Parameters mirror `operon_fig_patched`.

    Returns
    -------
    pandas.DataFrame
    """
    out = df.copy()
    out['ID'] = out[group_col] if group_col in out.columns else 'Unknown_Block'
    out['org_name'] = out[org_col] if org_col in out.columns else 'Unknown Organism'
    out = resolve_domain_labels(out, label_col)
    out = flag_query_rows(out)
    out['pid_order'] = pd.factorize(out['ID'])[0]
    return out.reset_index(drop=True)


def compute_label_width(df):
    """
    Character width every row's 3-line label should be padded to, so the
    query-id / block-id / organism-name boxes line up across rows.
    """
    max_id_len = df['ID'].astype(str).str.len().max()
    max_org_len = df['org_name'].astype(str).str.len().max()
    query_pids = df.loc[df['is_query'], 'pid'].astype(str)
    max_q_len = query_pids.str.len().max() if not query_pids.empty else 8
    return max(max_id_len, max_org_len, max_q_len)


# ---------------------------------------------------------------------------
# Row label
# ---------------------------------------------------------------------------

def pad_and_escape(text, width):
    """
    Right-pad `text` to `width` characters, HTML-escape it, then turn the
    padding spaces into non-breaking spaces (&nbsp;).

    Graphviz HTML-like labels collapse normal spaces, so padding only
    works if it survives as &nbsp;. Doing this with a monospace font
    (see `build_row_label_html`) is what makes the label boxes line up
    across rows.
    """
    return html.escape(str(text).ljust(width)).replace(" ", "&nbsp;")


def build_row_label_html(query_pid, block_id, org_name, width, font_size):
    """
    Build the 3-line Graphviz HTML-like label for one block's left-hand
    label box: bold query protein id, block id, italic organism name.
    """
    query_str = pad_and_escape(query_pid, width)
    block_str = pad_and_escape(block_id, width)
    org_str = pad_and_escape(org_name, width)
    return (
        '<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
        f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{font_size}">'
        f'<B>{query_str}</B></FONT></TD></TR>'
        f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{font_size}">'
        f'{block_str}</FONT></TD></TR>'
        f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas italic" POINT-SIZE="{font_size}">'
        f'{org_str}</FONT></TD></TR>'
        '</TABLE>>'
    )


# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------

def build_color_map(df, max_colors=5, ignore_domains=None, custom_colors=None):
    """
    Decide which fill color each domain/annotation value gets.

    Priority order:

    1. Anything explicitly listed in `custom_colors` keeps that exact
       color. These do not count against `max_colors`.
    2. Domains seen on a query gene get an automatic color next (these
       are usually the reason the figure exists).
    3. The remaining `max_colors` slots go to the most frequent of the
       remaining domains, by gene count.
    4. Everything else -- including anything in `ignore_domains`, the
       literal value 'unk'/blank/'-'/'?', and anything containing the
       word "hypothetical" -- is left white (no color).

    Parameters
    ----------
    df : pandas.DataFrame
        Must already have 'domain' and 'is_query' columns (see
        `prepare_dataframe`).
    max_colors : int
        Number of *automatically chosen* colors, on top of any fixed via
        `custom_colors`.
    ignore_domains : list[str] or None
        Domain values that never get an automatic color (case
        insensitive). Defaults to `DEFAULT_IGNORE_DOMAINS`.
    custom_colors : dict[str, str] or None
        Explicit {domain_value: color} overrides. Keys should match the
        values that appear in the 'domain' column -- i.e. typically the
        raw Pfam identifiers when `label_col='pfam'` (the default).
        Values can be any color Graphviz/seaborn understands, e.g.
        '#ff8800' or 'tomato'.
        Example: `custom_colors={'LysR_substrate': '#ff8800', 'PrpF': 'tomato'}`

    Returns
    -------
    dict[str, str]
        Mapping from domain value to color string.
    """
    if ignore_domains is None:
        ignore_domains = DEFAULT_IGNORE_DOMAINS
    custom_colors = dict(custom_colors) if custom_colors else {}

    ignore_lower = [d.lower() for d in ignore_domains] + ['unk', ' ', '-', '?']
    domain_str = df['domain'].astype(str)
    is_ignorable = (
        domain_str.str.lower().isin(ignore_lower)
        | domain_str.str.lower().str.contains('hypothetical', na=False)
    )
    already_colored = domain_str.isin(custom_colors)

    query_domains = (
        df.loc[df['is_query'] & ~is_ignorable & ~already_colored, 'domain']
        .unique()
        .tolist()
    )

    remaining_slots = max(0, max_colors - len(query_domains))
    freq_domains = (
        df.loc[~is_ignorable & ~already_colored & ~domain_str.isin(query_domains), 'domain']
        .value_counts()
        .head(remaining_slots)
        .index
        .tolist()
    )

    auto_domains = query_domains + freq_domains
    palette = sns.color_palette('pastel', len(auto_domains)).as_hex() if auto_domains else []

    color_map = dict(custom_colors)
    color_map.update(dict(zip(auto_domains, palette)))
    return color_map


# ---------------------------------------------------------------------------
# Orientation
# ---------------------------------------------------------------------------

def normalize_block_strand(block_df, normalize_orientation=True):
    """
    Optionally mirror a block so its query gene always points the same
    way (strand +1, drawn as a right-pointing arrow).

    Neighborhoods are usually pulled out with no regard for which strand
    the query happens to land on, which makes "upstream"/"downstream"
    mean different things from row to row. When `normalize_orientation`
    is True and this block's query is on the minus strand, the whole
    block is reversed -- gene order *and* every strand sign flip -- which
    is equivalent to flipping the picture so it reads in the same
    direction as every other block.

    If a block has more than one query row, the first one is used to
    decide whether to flip.

    Parameters
    ----------
    block_df : pandas.DataFrame
        Rows for a single block, already in genomic (left-to-right) order.
    normalize_orientation : bool

    Returns
    -------
    pandas.DataFrame
        `block_df` unchanged, or a reversed/strand-flipped copy.
    """
    if not normalize_orientation:
        return block_df

    query_rows = block_df[block_df['is_query']]
    if query_rows.empty or query_rows['strand'].iloc[0] != -1:
        return block_df

    flipped = block_df.iloc[::-1].copy()
    flipped['strand'] = -flipped['strand']
    return flipped


def collapsed_triangle_orientation(query_canonical_strand):
    """
    Graphviz `orientation` (degrees) that makes a `shape=triangle` node
    point in the direction *opposite* the query's canonical direction.

    A plain triangle points up by default. Empirically (verified by
    rendering test shapes with Graphviz 2.43): rotating it 90 degrees
    points the apex left, and 270 (equivalently -90) points it right.
    So if the query points right (strand +1, drawn as `rarrow`), an
    opposite-strand gene collapses to a left-pointing triangle, and
    vice versa.
    """
    return 90 if query_canonical_strand == 1 else 270


# ---------------------------------------------------------------------------
# Per-gene node styling
# ---------------------------------------------------------------------------

def gene_node_style(row, query_canonical_strand, color_map, highlight_query=True,
                     collapse_opposite_strand=False, font_size=10):
    """
    Decide the Graphviz node attributes (shape/color/label/size) for one
    gene.

    By default every gene is drawn at full size: an arrow (right for
    strand +1, left for strand -1, a plain box for anything else),
    filled per `color_map` and labeled with its domain/annotation. The
    query gene additionally gets a red, thicker outline when
    `highlight_query` is True.

    When `collapse_opposite_strand` is True, genes on the strand
    *opposite* the query's (`query_canonical_strand`) are drawn instead
    as small, unlabeled, grey triangles pointing away from the query's
    direction -- a lightweight "something is here, transcribed the
    other way" cue instead of giving them the same visual weight as
    genes that share the query's reading direction. The query gene
    itself is never collapsed.

    Parameters
    ----------
    row : pandas.Series
        One gene row; needs 'strand', 'domain' and 'is_query'.
    query_canonical_strand : int
        The strand (1 or -1) this block's query is drawn with.
    color_map : dict[str, str]
    highlight_query : bool
    collapse_opposite_strand : bool
    font_size : int

    Returns
    -------
    dict
        Keyword arguments for `AGraph.add_node`.
    """
    strand_val = row.get('strand', 1)
    is_target = bool(row['is_query'])

    opposite_strand = (
        collapse_opposite_strand
        and not is_target
        and strand_val != query_canonical_strand
    )

    if opposite_strand:
        return dict(
            label='',
            shape='triangle',
            orientation=collapsed_triangle_orientation(query_canonical_strand),
            style='filled',
            fixedsize='true',
            width='0.18',
            height='0.18',
            fillcolor='#cccccc',
            color='#888888',
            penwidth='1',
        )

    node_shape = 'rarrow' if strand_val == 1 else 'larrow' if strand_val == -1 else 'box'
    return dict(
        label=str(row['domain']),
        shape=node_shape,
        style='filled',
        fixedsize='false',  # text length dictates the box size naturally
        margin='0.1,0.05',
        height='0.4',
        color='red' if (highlight_query and is_target) else 'black',
        penwidth='3' if (highlight_query and is_target) else '1',
        fillcolor=color_map.get(row['domain'], '#ffffff'),
        fontsize=font_size,
        fontname='Consolas',
    )


# ---------------------------------------------------------------------------
# Graph assembly
# ---------------------------------------------------------------------------

def add_block_to_graph(graph, block_df, block_index, label_width, color_map,
                        highlight_query=True, collapse_opposite_strand=False,
                        font_size=10, left_pad=0, right_pad=0, spacer_width=0.6):
    """
    Add one full row (label box + every gene) to `graph`, wiring
    everything together with invisible same-rank edges so it is drawn as
    a single left-to-right row.

    `left_pad`/`right_pad` invisible spacer nodes are inserted before the
    first gene / after the last gene respectively, so that rows with
    fewer genes than the widest row still take up the same amount of
    horizontal space on each side of the query. This is what lets
    `operon_fig_patched(..., align_query_center=True)` line the query
    gene up in (approximately) the same column on every row -- "approximately"
    because real gene boxes have label-dependent widths, so this is a
    layout heuristic, not a pixel-exact guarantee. See `chain_align_nodes`
    for the second part of that trick.

    Parameters
    ----------
    graph : pygraphviz.AGraph
    block_df : pandas.DataFrame
        Rows for one block, in the left-to-right order to draw them in
        (already normalized/reversed by `normalize_block_strand` if
        that was requested).
    block_index : int
        Used to build unique node ids for this row.
    label_width : int
        From `compute_label_width`.
    color_map : dict[str, str]
    highlight_query, collapse_opposite_strand, font_size :
        Forwarded to `gene_node_style`.
    left_pad, right_pad : int
        Number of invisible spacer nodes to add on each side.
    spacer_width : float
        Width (inches) of each spacer node.

    Returns
    -------
    dict
        {'label_node': node id, 'query_node': node id or None,
         'gene_nodes': [node ids, left-to-right]}
    """
    query_rows = block_df[block_df['is_query']]
    query_pid = query_rows['pid'].iloc[0] if not query_rows.empty else 'No Query'
    query_canonical_strand = query_rows['strand'].iloc[0] if not query_rows.empty else 1

    label_node_id = f'label_{block_index}'
    label_html = build_row_label_html(
        query_pid=query_pid,
        block_id=block_df['ID'].iloc[0],
        org_name=block_df['org_name'].iloc[0],
        width=label_width,
        font_size=font_size,
    )
    graph.add_node(label_node_id, label=label_html, shape='none', margin=0.1)

    gene_node_ids = []
    query_node_id = None
    for row_position, (_, row) in enumerate(block_df.iterrows()):
        node_id = f'gene_{block_index}_{row_position}'
        style = gene_node_style(
            row,
            query_canonical_strand=query_canonical_strand,
            color_map=color_map,
            highlight_query=highlight_query,
            collapse_opposite_strand=collapse_opposite_strand,
            font_size=font_size,
        )
        graph.add_node(node_id, **style)
        gene_node_ids.append(node_id)
        if row['is_query']:
            query_node_id = node_id

    spacer_ids_left = [f'spacer_{block_index}_L{i}' for i in range(left_pad)]
    spacer_ids_right = [f'spacer_{block_index}_R{i}' for i in range(right_pad)]
    for spacer_id in spacer_ids_left + spacer_ids_right:
        graph.add_node(spacer_id, label='', shape='box', style='invis',
                        width=spacer_width, height=0.01)

    row_node_ids = [label_node_id] + spacer_ids_left + gene_node_ids + spacer_ids_right
    graph.add_subgraph(row_node_ids, rank='same')
    for a, b in zip(row_node_ids[:-1], row_node_ids[1:]):
        graph.add_edge(a, b, style='invis', penwidth=0)

    return {'label_node': label_node_id, 'query_node': query_node_id, 'gene_nodes': gene_node_ids}


def chain_align_nodes(graph, node_ids, weight=10000):
    """
    Pull a list of nodes -- one per row, in row order -- into the same
    visual column, by connecting consecutive nodes with a very high
    weight invisible edge.

    `dot` lays a graph out by minimizing total (edge weight x edge
    length); putting a large weight on an edge that is never actually
    drawn (`style='invis'`) is the standard trick to bias two nodes on
    different ranks towards the same horizontal position. The original
    code already used this for the row-label boxes; this function pulls
    that out so it can be reused for the query column too.

    Parameters
    ----------
    graph : pygraphviz.AGraph
    node_ids : list
        Skips silently over any `None` entries (e.g. a block with no
        query gene).
    weight : int
    """
    clean_ids = [n for n in node_ids if n is not None]
    for a, b in zip(clean_ids[:-1], clean_ids[1:]):
        graph.add_edge(a, b, style='invis', penwidth=0, weight=weight)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def operon_fig_patched(df, group_col='block_id', label_col='pfam', org_col='organism',
                        output_file='operon_fig_out.svg', max_colors=5,
                        highlight_query=True, font_size=10, ignore_domains=None,
                        custom_colors=None, normalize_orientation=True,
                        align_query_center=True, collapse_opposite_strand=False,
                        spacer_width=0.6):
    """
    Draw a gene-neighborhood ("operon") figure, one row per block, and
    write it to `output_file`.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-format table, one row per gene. See the module docstring
        for required columns.
    group_col : str
        Column that identifies which block/row a gene belongs to.
    label_col : str
        Column used for node labels/colors. 'pfam' (the default) uses
        the pfam -> aravind -> product fallback chain; see
        `resolve_domain_labels`.
    org_col : str
        Column with the organism name shown in the row label.
    output_file : str
        Path Graphviz writes the rendered figure to (extension controls
        format, e.g. '.svg', '.png', '.pdf').
    max_colors : int
        Max number of *automatically* assigned colors; see `build_color_map`.
    highlight_query : bool
        Outline the query gene in red.
    font_size : int
    ignore_domains : list[str] or None
        Domain values to never auto-color (see `build_color_map`).
    custom_colors : dict[str, str] or None
        Explicit {domain_value: color} overrides -- keys should match
        values in the 'domain' column (typically raw pfam ids when
        `label_col='pfam'`). See `build_color_map`.
    normalize_orientation : bool, default True
        If True, blocks whose query is on the minus strand are mirrored
        so every query is drawn pointing the same way (strand +1).
        Set to False to draw every block in its original orientation
        (the old, pre-this-change behavior).
    align_query_center : bool, default True
        If True, pad each row with invisible spacer nodes (and an extra
        alignment edge) so the query gene falls in roughly the same
        column on every row, instead of each row simply starting flush
        left. See `add_block_to_graph` and `chain_align_nodes` for the
        mechanics and its limits.
    collapse_opposite_strand : bool, default False
        If True, genes on the opposite strand from the query are drawn
        as small unlabeled grey triangles pointing away from the
        query's direction, instead of full domain-labeled arrows. Useful
        for decluttering when those genes are not of interest. The query
        gene itself is never collapsed.
    spacer_width : float, default 0.6
        Width (inches) of the invisible spacer nodes used for
        `align_query_center`. Tune this if real gene boxes in your
        figure are consistently much narrower/wider than this.

    Returns
    -------
    pandas.DataFrame
        The working copy of `df` actually used to build the figure
        (with 'ID', 'org_name', 'domain', 'is_query', 'pid_order' added),
        mainly useful for debugging.
    """
    working = prepare_dataframe(df, group_col=group_col, org_col=org_col, label_col=label_col)
    label_width = compute_label_width(working)
    color_map = build_color_map(
        working, max_colors=max_colors, ignore_domains=ignore_domains, custom_colors=custom_colors
    )

    blocks = [
        normalize_block_strand(block_df, normalize_orientation=normalize_orientation)
        for _, block_df in working.groupby('pid_order', sort=True)
    ]

    # First pass: how far left/right of the query does each block extend?
    # Needed up front so every row can be padded to the same width before
    # any nodes are added.
    left_counts, right_counts = [], []
    for block_df in blocks:
        query_positions = np.flatnonzero(block_df['is_query'].to_numpy())
        q_pos = query_positions[0] if len(query_positions) else 0
        left_counts.append(q_pos)
        right_counts.append(len(block_df) - 1 - q_pos)
    max_left = max(left_counts) if (align_query_center and left_counts) else 0
    max_right = max(right_counts) if (align_query_center and right_counts) else 0

    graph = pgv.AGraph(directed=True)
    graph.graph_attr.update(nodesep=0.05, ranksep=0.15)

    blocks_info = []
    for block_index, (block_df, n_left, n_right) in enumerate(zip(blocks, left_counts, right_counts)):
        left_pad = (max_left - n_left) if align_query_center else 0
        right_pad = (max_right - n_right) if align_query_center else 0
        info = add_block_to_graph(
            graph, block_df, block_index, label_width, color_map,
            highlight_query=highlight_query,
            collapse_opposite_strand=collapse_opposite_strand,
            font_size=font_size,
            left_pad=left_pad,
            right_pad=right_pad,
            spacer_width=spacer_width,
        )
        blocks_info.append(info)

    # Align the label boxes into one column (original behavior)...
    chain_align_nodes(graph, [info['label_node'] for info in blocks_info])
    # ...and, optionally, the query genes into their own column.
    if align_query_center:
        chain_align_nodes(graph, [info['query_node'] for info in blocks_info])

    graph.draw(output_file, prog='dot')
    return working


if __name__ == '__main__':
    # Small smoke test using a couple of synthetic blocks: one query on
    # the plus strand, one on the minus strand (to exercise orientation
    # normalization), each with a mix of same- and opposite-strand
    # neighbors (to exercise collapsing).
    demo = pd.DataFrame([
        # block A: query on + strand
        dict(block_id='blockA', organism='Org A', pid='A_1', strand=1, query=0, pfam='DomX'),
        dict(block_id='blockA', organism='Org A', pid='A_2', strand=-1, query=0, pfam='DomY'),
        dict(block_id='blockA', organism='Org A', pid='A_3', strand=1, query=1, pfam='QueryDom'),
        dict(block_id='blockA', organism='Org A', pid='A_4', strand=1, query=0, pfam='DomZ'),
        # block B: query on - strand (should flip if normalize_orientation=True)
        dict(block_id='blockB', organism='Org B', pid='B_1', strand=-1, query=0, pfam='DomZ'),
        dict(block_id='blockB', organism='Org B', pid='B_2', strand=-1, query=1, pfam='QueryDom'),
        dict(block_id='blockB', organism='Org B', pid='B_3', strand=1, query=0, pfam='DomY'),
    ])
    operon_fig_patched(
        demo,
        output_file='demo_operon_fig.svg',
        custom_colors={'QueryDom': '#ff8800'},
        normalize_orientation=True,
        align_query_center=True,
        collapse_opposite_strand=True,
    )