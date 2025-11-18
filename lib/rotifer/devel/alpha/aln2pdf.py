from re import sub
from reportlab.lib.colors import lightblue, black
import pdfplumber
from reportlab.pdfgen import canvas
from reportlab.lib.colors import black, lightblue, Color
from PyPDF2 import PdfReader, PdfWriter
from weasyprint.pdf.stream import subset
pastel_orange = Color(1.0, 0.8, 0.6) 

def draw_arrow_from_to(
    c,
    x_start,
    x_end,
    y,
    shaft_height=6,
    arrow_proportion=0.2,
    arrow_width=None,
    color=pastel_orange,
    edge_color=black,
    edge_width=1,
    show_edge=True
):
    """
    Draws a horizontal arrow from x_start to x_end, vertically centered on y.
    - Arrowhead size is based on proportion of total length
    - Allows stroke (outline) control and default pastel orange fill
    """
    total_length = x_end - x_start
    print (total_length)
    if total_length < 50:
        arrow_proportion = 0.5 - (0.4 / (50 - 10)) * (total_length - 10)
    else:
        arrow_proportion = 0.1

    arrow_size = total_length * arrow_proportion
    shaft_end_x = x_end - arrow_size
    shaft_half = shaft_height / 2
    arrow_width = arrow_width or shaft_height * 2
    arrow_half_width = arrow_width / 2
    stroke_flag = 1 if show_edge else 0

    # Set colors and stroke style
    c.setFillColor(color)
    c.setStrokeColor(edge_color)
    c.setLineWidth(edge_width)

    # Draw shaft
    c.rect(x_start, y - shaft_half, shaft_end_x - x_start, shaft_height, fill=1, stroke=stroke_flag)

    # Draw arrowhead
    tip = (x_end, y)
    base_top = (shaft_end_x, y + arrow_half_width)
    base_bottom = (shaft_end_x, y - arrow_half_width)

    path = c.beginPath()
    path.moveTo(*tip)
    path.lineTo(*base_top)
    path.lineTo(*base_bottom)
    path.close()
    c.drawPath(path, fill=1, stroke=stroke_flag)


def draw_cylinder_from_to(
    c,
    x_start,
    x_end,
    y,
    height=4,
    color=lightblue,
    edge_color=black,
    edge_width=0.1,
    show_edge=True,
    label=None
):
    """
    Draws a horizontal 3D-style cylinder (helix) from x_start to x_end, centered vertically on y.
    Allows full control over the stroke (outline).
    """
    width = x_end - x_start
    r = height / 2
    y_bottom = y - r

    # Setup colors and stroke
    c.setFillColor(color)
    c.setStrokeColor(edge_color)
    c.setLineWidth(edge_width)

    stroke_flag = 1 if show_edge else 0

    # Draw main body
    c.roundRect(x_start, y_bottom, width, height, radius=r, fill=1, stroke=stroke_flag)
    c.setFillColor('#d48c82')

    # Rear cap (right)
    c.ellipse(x_end - r * 2, y_bottom, x_end, y_bottom + height, fill=1, stroke=stroke_flag)

    # Front cap (left)
    
    #c.ellipse(x_start, y_bottom, x_start + r * 2, y_bottom + height, fill=1, stroke=stroke_flag)

    # Optional label
    if label:
        c.setFillColor(black)
        c.setFont("Helvetica", 6)
        c.drawCentredString((x_start + x_end) / 2, y_bottom + height + 3, label)

def add_SS_diagrams(pdf_in, pdf_out):
    from func import draw_cylinder_from_to, draw_arrow_from_to
    with pdfplumber.open(pdf_in) as pdf:
        page = pdf.pages[0]
        width, height = page.width, page.height
        words = page.extract_words()
        l =[]
        for x in words:
            if x['text'] =='SS':
                for y in words:
                    if y['top'] == x['top']:
                        if y['text'] != x['text']:
                            l.append(y)




    # === STEP 2: Create overlay with ReportLab ===
    overlay_path = "overlay.pdf"
    c = canvas.Canvas(overlay_path, pagesize=(width, height))
    for x in l:
        # Draw centered arrow
        if x['text'].startswith('E'):
            draw_arrow_from_to(c, x_start=x['x0'], x_end=x['x1'], y=height - x['top'],arrow_proportion=0.10, shaft_height=2, arrow_width=4, edge_width=0.1)
        else:
        # Draw two cylinders (helices) beside the arrow
            draw_cylinder_from_to(c,x_start= x['x0'], x_end=x['x1'],y=height - x['top'], height=4, color=lightblue, edge_width=0.1)


    c.save()

    # === STEP 3: Merge overlay into original PDF ===
    reader = PdfReader(pdf_in)
    overlay = PdfReader(overlay_path)
    writer = PdfWriter()

    base_page = reader.pages[0]
    overlay_page = overlay.pages[0]
    base_page.merge_page(overlay_page)
    writer.add_page(base_page)

    with open(pdf_out, "wb") as f:
        writer.write(f)

    print(f"✅ Done! Annotated PDF saved to: {pdf_out}")



def veremos(in_aln_r,
            consensus=False,
            annotations=False,
            ss=False,
            font_size=6,
            output='figure.pdf',
            landscape=False,
            aln_length=70,
            ss_as_polygon=False
             ):
    import numpy as np
    aln_r = in_aln_r.copy()
    from rotifer.core.functions import chunks
    from rotifer.devel.alpha.aln2pdf import html_highlight_aln, html_highlight_consensus, make_dummy_df
    from weasyprint import HTML
    import sys
    if not consensus:
        consensus = aln_r[aln_r.index.str.contains("Consensus", case=False)].index[0]
    if not ss:
        try:
            ss = aln_r[aln_r.index.str.startswith("@")].index.str.lstrip('@').tolist()
            aln_r.index = aln_r.index.str.lstrip('@')
        except:
            ss = False 
    if not annotations:
        try:
            annotations = aln_r[aln_r.index.str.startswith("#")].index.str.lstrip('#').tolist()
            aln_r.index = aln_r.index.str.lstrip('#')
        except:
            annotations = False 

    if landscape:
        landscape='landscape'
    else:
        landscape =""
    if ss_as_polygon:
        aln_r.loc[ss] = aln_r.loc[ss].replace('-',' ')
    headers = {
        'selector': 'th:not(.index_name)',
        'props': f'''font-size: {font_size}px;
        text-align: left;
        font-family:"Courier New";
        color:black;'''
    }
    slices = chunks(aln_r.columns, aln_length)
    # using slices to create a dummy df to later calculate the size when chunked 
    slices = list(slices)
    to_merge = []

    for x in slices:
        sliced = aln_r.loc[:,x].copy()
        slice_consensus = ([consensus], sliced.columns)
        if annotations:
            slice_annotaion = (annotations, sliced.columns)
            con_ann = annotations + [consensus]
        else:
            con_ann = [consensus]
        if ss:
            slice_ss = (ss, sliced.columns)
            con_ann = con_ann + ss
            

        slice_sequences = sliced.loc[~sliced.index.isin(con_ann)].index.tolist() 
        slice_sequences = slice_sequences + [consensus]
        slice_sequences = (slice_sequences, sliced.columns)
        if sys.version_info.minor > 8:
            df_style = sliced.loc[:,x].style.set_properties(**{
                'font-size': f'{font_size}px',
                'font-family':"'DejaVu Sans Mono', monospace",
                "text-align": "center"}
            ).apply(html_highlight_aln, axis=0, subset=slice_sequences).hide(axis='columns').apply(
                html_highlight_consensus, subset=slice_consensus
            ).set_table_styles(
                [headers]
            )
            if ss_as_polygon:
                df_style = df_style.apply(html_white_text, subset=slice_ss)

        else:
            df_style = sliced.loc[:,x].style.set_properties(**{
                'font-size': f'{font_size}px',
                'font-family':"'DejaVu Sans Mono', monospace",
                "text-align": "center"}
            ).apply(html_highlight_aln, axis=0, subset=slice_sequences).hide_columns().apply(
                html_highlight_consensus, subset=slice_consensus
            ).set_table_styles(
                [headers]
            )
            if ss_as_polygon:
                df_style = df_style.apply(html_white_text, subset=slice_ss)
                
        df_style = df_style.to_html(table_attributes='cellspacing=0, cellpadding=0')
        to_merge.append(df_style)
    ht = '<br>'.join(to_merge)        
    ##### Estimating the size to create a canvas using a dummmy df:

    dummydf = make_dummy_df(aln_r.shape[0] * len(slices), aln_length)
    dummydf.index = np.tile(aln_r.index, len(slices))
    w,h = estimate_canvas_size(dummydf, font_size=font_size, baseline_char_width=1.1, baseline_row_height=5.9)
    html_full = f"""
    <html>
    <head>
      <meta charset="utf-8">
      <style>
      @page {{
      size: {w}pt {h}pt;
      margin: 1cm;
        }}

        table, th, td {{
          border: none !important;
          border-collapse: collapse;
        }}
      </style>
    </head>
    <body>
      {ht}
    </body>
    </html>
    """
    HTML(string=html_full).write_pdf(f"{output}")



def polish_2_residues_df(polish_file):
    import re
    import pandas as pd
    from collections import Counter

    def drop_leading_empty_columns(df):
        while True:
            df.columns = range(df.shape[1])  # Reindex columns as 0, 1, 2, ...
            if df.empty or df.shape[1] == 0:
                break
            first_col = df.columns[0]
            # Check if entire first column is empty (NaN or blank string)
            if df[first_col].apply(lambda x: pd.isna(x) or (isinstance(x, str) and x.strip() == '')).all():
                df = df.drop(columns=first_col)
            else:
                break
        return df

    def insert_tabs(line, positions):
        chars = list(line)

        # Shift mode_second to AFTER the space
        adjusted_positions = []
        for i, pos in enumerate(positions):
            if pos is None:
                continue
            # If it's the second position (index 1), insert after the space
            adjusted_pos = pos + 1 if i == 1 else pos
            if 0 <= adjusted_pos <= len(line):
                adjusted_positions.append(adjusted_pos)

        # Insert tabs in reverse order to avoid shifting issues
        for pos in sorted(set(adjusted_positions), reverse=True):
            chars.insert(pos, '\t')

        return ''.join(chars)

    # Step 1: Read lines
    with open(polish_file, "r") as f:
        lines = [line.rstrip('\n') for line in f]

    # Step 2: Extract space-after-word positions per line

    results = [[m.start() for m in re.finditer(r'(?<=\S)\s', line)] for line in lines]

    #results = [[m.start(2) for m in re.finditer(r'(\w+)(\s+ )', line)] for line in lines]

    # Step 3: Compute stats across all lines
    first_elements = [lst[0] for lst in results if len(lst) >= 1]
    second_elements = [lst[1] for lst in results if len(lst) >= 2]
    last_elements = [lst[-1] for lst in results if lst]

    max_first = max(first_elements) if first_elements else None
    mode_second = Counter(second_elements).most_common(1)
    mode_second = mode_second[0][0] if mode_second else None
    mode_last = Counter(last_elements).most_common(1)
    mode_last = mode_last[0][0] if mode_last else None

    # Step 4: Function to insert tabs at multiple positions

    # Step 5: Apply to each line
    new_lines = [insert_tabs(line, [max_first, mode_second, mode_last]) for line in lines]
    tsv_string = '\n'.join(new_lines)
    # Step 6: Save result in pd.DataFrame
    from io import StringIO
    df = pd.read_csv(StringIO(tsv_string), sep='\t', header=None)
    rdf = df.set_index([0])[2].astype(str).str.split("", expand=True)
    rdf = drop_leading_empty_columns(rdf)
    rdf.insert(0,'start',df[1].astype(str).str.strip().values)
    rdf['end'] = df[3].astype(str).str.strip().fillna('').values
    rdf.columns = list((range(1,rdf.shape[1]+1)))
    rdf.index = rdf.index.str.strip()
    rdf = rdf.replace('', ' ')
    rdf.index.name =None

    return rdf

def html_highlight_aln(s):
    from rotifer.core.functions import loadConfig
    from rotifer.core  import config as CoreConfig
    import numpy as np
    import pandas as pd 
    cd = loadConfig(
            ':colors.html_aa_colors',
            system_path=CoreConfig['baseDataDirectory'])
    ### getting the consensus value to map the colors filling na with "  " to color white 
    d = cd[s.fillna('_').iloc[-1]]
    #d = aa_groups_colors[s.fillna('  ').iloc[-1]]
    return np.where(
        s == '  ',
        'color:"";background-color:""',
        np.where(
            s == s.iloc[-1],
            f'color:{d["fcolor"]};background-color:{d["color"]};font-weight: bold;text-align: center',
            np.where(pd.to_numeric(s, errors="coerce").notna(),
                     'color:green;background-color:',np.where(
                        s.isin(d['residues']),
                            f'color:{d["fcolor"]};background-color:{d["color"]};font-weight: bold;text-align: center',
                            f'color:#74718c;background-color:;font-weight: bold;text-align: center'))))

def html_highlight_consensus(s):
    from rotifer.core.functions import loadConfig
    from rotifer.core  import config as CoreConfig
    from rotifer.devel.alpha import gian_func as gf
    import numpy as np
    cd = loadConfig(
            ':colors.html_aa_colors',
            system_path=CoreConfig['baseDataDirectory'])
    d = cd[s.fillna('_').iloc[-1]]
    """TODO: Docstring for highlight_consensus.

    :arg1: TODO
    :returns: TODO

    """
    return np.where(
        s.isin(cd["ALL"]["residues"]),
        f'color:{d["fcolor"]};background-color:{d["color"]};font-weight: bold;text-align: center',
        f'color:{d["fcolor"]};background-color:{d["color"]};font-weight: bold;text-align: center',
        )
def html_white_text(s):
    import pandas as pd
    return pd.Series(["color: white; font-family: Courier New"] * len(s), index=s.index)

def clean_and_rename_columns(df):
    import pandas as pd
    # Drop columns where all values are NaN, "" or " "
    non_empty = df.applymap(lambda x: not (pd.isna(x) or str(x).strip() == "")).any()
    df_cleaned = df.loc[:, non_empty]

    # Rename columns to 1, 2, 3, ...
    df_cleaned.columns = range(1, df_cleaned.shape[1] + 1)

    return df_cleaned

def add_SS_diagrams_pages(pdf_in, pdf_out):
    from rotifer.core.functions import loadConfig
    from rotifer.core  import config as CoreConfig
    color_dict = loadConfig(
            ':colors.pdf_SS_colors',
            system_path=CoreConfig['baseDataDirectory'])
    from rotifer.devel.alpha.aln2pdf import draw_cylinder_from_to, draw_arrow_from_to
    from io import BytesIO
    reader = PdfReader(pdf_in) # added later on
    output = PdfWriter()
    with pdfplumber.open(pdf_in) as pdf:
        for i,page in enumerate(pdf.pages):
            width, height = page.width, page.height

            words = page.extract_words()
            l =[]
            for x in words:
                if x['text'].startswith('Secondary'):
                    for y in words:
                        if y['top'] == x['top']:
                            if y['text'] != x['text']:
                                l.append(y)

        # === STEP 2: Create overlay with ReportLab ===
            packet = BytesIO()
            c = canvas.Canvas(packet, pagesize=(width, height))
            for x in l:
                # Draw centered arrow
                if x['text'].startswith('E'):
                    draw_arrow_from_to(c,
                                       x_start=x['x0'],
                                       x_end=x['x1'],
                                       y=height - x['top'],
                                       arrow_proportion=0.10,
                                       shaft_height=4,
                                       arrow_width=8,
                                       edge_width=0.1,
                                       color=color_dict['strand']['color'],
                                       edge_color=color_dict['strand']['line'])
                else:
                # Draw two cylinders (helices) beside the arrow
                    draw_cylinder_from_to(c,
                                          x_start= x['x0'],
                                          x_end=x['x1'],
                                          y=height - x['top'],
                                          height=4,
                                          color=color_dict['helix']['color'],
                                          edge_color=color_dict['helix']['line'],
                                          edge_width=0.1)
            c.save()
            packet.seek(0)

            # === STEP 3: Merge overlay into original PDF ===
            overlay_pdf = PdfReader(packet)
            #original_page = PdfReader(pdf_in).pages[i] Trying to fix issues 
            original_page = reader.pages[i]
            original_page.merge_page(overlay_pdf.pages[0])
            output.add_page(original_page)


    with open(pdf_out, "wb") as f:
        output.write(f)

    print(f"✅ Done! Annotated PDF saved to: {pdf_out}")


def estimate_canvas_sizei2(df, font_size, baseline_font_size=4.5, baseline_char_width=1.5, baseline_row_height=5.2):
    """
    Estimate canvas size in pixels for a given font size, scaled from a known baseline (4.5pt).

    Parameters:
    - df: pandas DataFrame
    - font_size: the font size you plan to use
    - baseline_font_size: reference font size (4.5 pt)
    - baseline_char_width: reference character width in pixels at baseline (1.5 px)
    - baseline_row_height: reference row height in pixels at baseline (5.2 px)

    Returns:
    - (total_width_px, total_height_px)
    """
    scale = font_size / baseline_font_size
    char_width_px = baseline_char_width * scale
    row_height_px = baseline_row_height * scale

    max_chars_per_col = [
        max(len(str(cell)) for cell in [col] + df[col].tolist())
        for col in df.columns
    ]

    total_width = sum(n * char_width_px for n in max_chars_per_col)
    total_height = row_height_px * (len(df) + 1)

    return int(total_width), int(total_height)

def estimate_canvas_size(df, font_size,
                           baseline_font_size=4.5,
                           baseline_char_width=1.2,
                           baseline_row_height=5.9):
    """
    Estimate canvas size in pixels for a given font size, scaled from a known baseline (4.5pt).
    Now includes index width in total width estimation.

    Parameters:
    - df: pandas DataFrame
    - font_size_pt: the font size you plan to use
    - baseline_font_size: reference font size (4.5 pt)
    - baseline_char_width: reference character width in pixels at baseline
    - baseline_row_height: reference row height in pixels at baseline

    Returns:
    - (total_width_px, total_height_px)
    """
    scale = font_size / baseline_font_size
    char_width_px = baseline_char_width * scale
    row_height_px = baseline_row_height * scale

    # Estimate max character count per column (including header)
    max_chars_per_col = [
        max(len(str(cell)) for cell in [col] + df[col].tolist())
        for col in df.columns
    ]

    # Estimate index width
    index_width = max(len(str(i)) for i in df.index)

    # Total width: index + all columns
    total_width = (index_width + sum(max_chars_per_col)) * char_width_px

    # Total height: one row per data entry, plus header
    total_height = row_height_px * (len(df) + 1)

    return int(total_width), int(total_height)

def make_dummy_df(n_rows, n_cols, fill="a"):
    import pandas as pd
    """
    Create a dummy DataFrame with the given shape, filled with a constant string.

    Parameters:
    - n_rows: number of rows
    - n_cols: number of columns
    - fill: string or value to fill each cell with

    Returns:
    - pandas DataFrame
    """
    return pd.DataFrame(
        [[fill] * n_cols for _ in range(n_rows)],
        columns=[f"col{i}" for i in range(n_cols)]
    )

