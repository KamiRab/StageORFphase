#!/usr/bin/env python3

""" Module handling all funtions related to plots with bokeh library. """

# Third
from pandas import Series, DataFrame
from numpy import nan, pi, mean
from bokeh.plotting import figure
from bokeh.models.glyphs import Quad
from bokeh.models import ColumnDataSource, Range1d, FuncTickFormatter
#from bokeh.charts import Bar, output_file, show
#from bokeh.charts.attributes import cat, color
#from bokeh.charts.operations import blend
from bokeh.layouts import gridplot
from bokeh.palettes import Spectral11, Category20, PuOr4, Category20c
from utils import poscov_to_dataframe
from utils import poscov_frames_signal
from utils import init_poscov
from utils import add_kmers_columns

# Local
from constants import *

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"


# Constants exclusive to this module
# Constants
PALETTE = Category20c[15]
PHASES_COLORS = {
    0:Category20[7][0], #blue
    1:Category20[7][-3],#green
    2:Category20[7][-1] #red
}
STRANDS_COLORS = {
    "+":"black",
    "-":"purple"
}
NT_COLORS = {
    "A":"green",
    "T":"red",
    "C":"blue",
    "G":"orange"
}
CODONS_COLORS = {
    "TAG":"red",
    "TAA":"red",
    "TGA":"red",
    "ATG":"green",
    "other":"#d9d9d9"
}
REG_COLORS = {
    CDS:PuOr4[2],
    UTR5:PuOr4[1],
    UTR3:PuOr4[1],
}
FRAMES_PAL = {
    0:PALETTE,
    1:PALETTE,
    2:PALETTE
}


def envs_plot(meta_start, meta_stop, bamfile, title="Envs",
              only_total=False, palette=PALETTE,
              colors=["blue", "red"]):
    """ Return a figure with both meta start and stop as VBar. """
    df_start = DataFrame(meta_start).transpose()
    df_stop = DataFrame(meta_stop).transpose()
    start_fig = vbar_kmers_figure(df_start, only_total=only_total,
                                  palette=palette, color=colors[0],
                                  title="Start")
    start_fig.xaxis.axis_label = "Relative positions to Start codon"
    stop_fig = vbar_kmers_figure(df_stop, only_total=only_total,
                                 palette=palette, color=colors[1],
                                 title="Stop")
    stop_fig.xaxis.axis_label = "Relative positions to Stop codon"
    # Link plots
    max_start = max([meta_start[i]["Total"] for i in sorted(meta_start)])
    max_stop = max([meta_stop[i]["Total"] for i in sorted(meta_stop)])
    max_y = max(max_start, max_stop)
    start_fig.y_range = Range1d(0, max_y)
    stop_fig.y_range = start_fig.y_range

    grid = gridplot([start_fig, stop_fig], ncols=1, responsive=True,
                    plot_width=600, plot_height=250)
    return grid

def plot_expression(gene_counts, types=[], title="VBar"):
    """ Plot VBar for each gene value in gene_counts.
    Types is a list where type of each gene is defined.
    NEED TE REIMPLEMENT TO TAKE IT FROM GFF TAG...
    """
    PALETTE = Category20[20]
    df = DataFrame({"counts":gene_counts})
    df["types"] = types

    colors = [*[PALETTE[0]]*6,
              *[PALETTE[2]]*5,
              PALETTE[4],
              *[PALETTE[6]]*3,
              PALETTE[8],
              PALETTE[10],
              *[PALETTE[12]]*9,
              *[PALETTE[14]]*7
              ]

    s = ColumnDataSource(df)
    bar_opts = dict(width=0.7, alpha=1)
    p = figure(title=title, plot_width=800, plot_height=350, responsive=True,
               x_range=list(df.index))
    p.vbar(bottom=0, top="counts", x=df.index,
           color=colors, legend="types", **bar_opts, source=s)
    p.yaxis.axis_label = "RPFs counts normalised by qPCR abundances"
    p.xaxis.major_label_orientation = pi/4
    return p

def vbar_kmers_figure(df, only_total=False, palette=PALETTE,
                      color="blue", title="VBars"):
    """ simple vbar figure with kmers stacked, or total bars. """
    bar_opts = dict(width=0.7, alpha=1)
    p = figure(title=title, plot_width=800, plot_height=350, responsive=True,
               x_range=Range1d(min(df.index), max(df.index)))

    if not only_total:
        sorted_kmers = sorted([k for k in df.columns if k != TOTAL_KEYNAME])
        top_cumul = Series()
        for (idx, kmer) in enumerate(sorted_kmers):
            kmer_color = palette[idx]
            # In this function, total is ignored to avoid 2x coverage
            if idx == 0:
                bottom = 0
                top = df[kmer]
            else:
                bottom = top_cumul
                top = bottom + df[kmer]
            # Adjut top dividing by the sizes of features
            p.vbar(bottom=bottom, top=top, x=sorted(df.index),
                   color=kmer_color, legend=str(kmer), **bar_opts)
            if top_cumul.empty:
                top_cumul = df[kmer]
            else:
                top_cumul = bottom + df[kmer]
    else:
        p.vbar(bottom=0, top=df[TOTAL_KEYNAME], x=sorted(df.index),
               color=color, **bar_opts)
    return p




def frames_plot(frames_counts, bamfile, title="Frames",
                only_total=False, palette=PALETTE,
                colors=["red", "blue", "green"]):
    """ Plot the Kmer repartition in each region of the transcripts. """
    bar_opts = dict(width=0.7, alpha=1)
    df = DataFrame(frames_counts).transpose()
    p = figure(title=title, plot_width=800, plot_height=350, responsive=True,
               x_range=["0", "1", "2"])

    if not only_total:
        sorted_kmers = sorted([k for k in df.columns if k != TOTAL_KEYNAME])
        top_cumul = Series()
        for (idx, kmer) in enumerate(sorted_kmers):
            kmer_color = palette[idx]
            # In this function, total is ignored to avoid 2x coverage
            if idx == 0:
                bottom = 0
                top = df[kmer]
            else:
                bottom = top_cumul
                top = bottom + df[kmer]
            # Adjut top dividing by the sizes of features
            p.vbar(bottom=bottom, top=top, x=[str(f) for f in df.index],
                   color=kmer_color, legend=str(kmer), **bar_opts)
            if top_cumul.empty:
                top_cumul = df[kmer]
            else:
                top_cumul = bottom + df[kmer]
    else:
        p.vbar(bottom=0, top=df[TOTAL_KEYNAME], x=[str(f) for f in df.index],
               color=colors, **bar_opts)
    return p

def regions_stats_plot(regs_stats, bamfile, title="Regs_stats",
                       only_total=False, palette=PALETTE,
                       colors=["black", "gray", "gray"]):
    """ Plot the Kmer repartition in each region of the transcripts. """
    bar_opts = dict(width=0.7, alpha=1)
    df = DataFrame(regs_stats)
    p = figure(title=title, plot_width=800, plot_height=350, responsive=True,
               x_range=ORDERED_REGIONS)
    kmers = add_kmers_columns(df, colname="counts")
    total_percentages = df[TOTAL_KEYNAME].sum()

    if not only_total:
        sorted_kmers = sorted([k for k in kmers if k != TOTAL_KEYNAME])
        top_cumul = Series()
        for (idx, kmer) in enumerate(sorted_kmers):
            str_kmer = str(kmer)
            kmer_percentages = (df[str_kmer] * 100) / total_percentages
            kmer_color = palette[idx]
            if idx == 0:
                bottom = 0
                top = kmer_percentages
            else:
                bottom = top_cumul
                top = bottom + kmer_percentages
            # Adjut top dividing by the sizes of features
            p.vbar(bottom=bottom, top=top, x=df.index,
                   color=kmer_color, legend=str_kmer, **bar_opts)
            if top_cumul.empty:
                top_cumul = kmer_percentages
            else:
                top_cumul = bottom + kmer_percentages
    else:
        percentages = (df[TOTAL_KEYNAME] * 100) / total_percentages
        p.vbar(bottom=0, top=percentages, x=df.index,
               color=colors, **bar_opts)

    p.yaxis.axis_label = "Percentages of total transcripts mapped reads"
    return p


def figure_splitted_signal(pos_cov,
                           transcript,
                           glyphs=["annot"],
                           sequence="",
                           title="Split frames",
                           only_total=False,
                           color=["blue", "green", "red"]):
    """ Return a grid with the signal of the 3 frames splitted. """
    ref_frame = transcript.CDS.start % 3
    frames_signals = {}
    for i in range(3):
        frames_signals[i] = init_poscov()
    # populate the frames_signals + return the maximum of the 3 frames
    frames_max_y = poscov_frames_signal(frames_signals, ref_frame, pos_cov)

    plots = []
    for frame in frames_signals:
        frame_poscov = frames_signals[frame]
        frame_palette = FRAMES_PAL[frame]
        frame_color = color[frame]
        if frame == 0:
            title = "frame %s (%s)" % (frame, transcript.gene_name)
        else:
            title = "frame %s" % (frame)
        frame_fig = figure_from_poscov(frame_poscov, transcript,
                                       glyphs=glyphs, sequence=sequence,
                                       title=title, palette=frame_palette,
                                       only_total=only_total,
                                       color=frame_color, max_y=frames_max_y)
        plots.append(frame_fig)

    # Adjut Y range
    plots[0].y_range.end = frames_max_y
    # Link all plots to the first frame
    plots[1].x_range = plots[0].x_range
    plots[2].x_range = plots[0].x_range
    plots[1].y_range = plots[0].y_range
    plots[2].y_range = plots[0].y_range
    # Only middle frame with Y axis name
    plots[0].yaxis.axis_label = ""
    plots[2].yaxis.axis_label = ""

    grid = gridplot(plots, ncols=1, responsive=True,
                    plot_width=600, plot_height=250)
    return grid


def figure_from_poscov(pos_cov, transcript, glyphs=["annot"], sequence="",
                       title="Figure from pos_cov", palette=Category20c[15],
                       only_total=False, color="black", max_y=None):
    """ Return a figure from the pos_cov dictionary. """
    pos_cov_df, found_kmers = poscov_to_dataframe(pos_cov)
    if only_total:
        tr_figure = total_vbar_figure(pos_cov_df, title=title, color=color)
    else:
        tr_figure = kmers_vbar_figure(pos_cov_df, found_kmers, title=title,
                                      palette=palette)
    if not max_y:
        max_y = max(pos_cov_df.coverages)
    else:
        pass
    add_annotations(
        tr_figure,
        transcript,
        max(pos_cov_df.index),
        max_y,
        glyphs=glyphs,
        sequence=sequence
        )
    return tr_figure

def total_vbar_figure(pos_cov_df, title="Total VBar", color="black"):
    """ Generate a figure with the total of read per position. """
    #print(pos_cov_df)
    bar_opts = dict(width=0.6, alpha=0.8)
    p = figure(title=title, plot_width=800, plot_height=350, responsive=True)
    p.vbar(bottom=0, top=pos_cov_df[TOTAL_KEYNAME].replace(0, nan),
           x=pos_cov_df.index, color=color, **bar_opts)
    p.xaxis.major_label_orientation = pi/4
    p.xaxis.formatter = FuncTickFormatter(code="""
    return %s[tick];
    """ % (list(pos_cov_df["positions"])))
    return p

def kmers_vbar_figure(pos_cov_df, kmers, title="VBar stacked kmers",
                      palette=Category20c[15]):
    """ Generate a figure with vbar for all kmers stacked. """
    bar_opts = dict(width=0.7, alpha=1)
    p = figure(title=title, plot_width=800, plot_height=350, responsive=True)
    sorted_kmers = sorted([k for k in kmers if k != TOTAL_KEYNAME])
    top_cumul = Series()
    for (idx, kmer) in enumerate(sorted_kmers):
        kmer_color = palette[idx]
        # In this function, total is ignored to avoid 2x coverage
        str_kmer = str(kmer)
        if idx == 0:
            bottom = 0
            top = pos_cov_df[str_kmer]
        else:
            bottom = top_cumul
            top = bottom + pos_cov_df[str_kmer]
        p.vbar(bottom=bottom, top=top, x=pos_cov_df.index,
               color=kmer_color, legend=str_kmer, **bar_opts)
        if top_cumul.empty:
            top_cumul = pos_cov_df[str_kmer]
        else:
            top_cumul = bottom + pos_cov_df[str_kmer]
    p.xaxis.major_label_orientation = pi/4
    p.xaxis.formatter = FuncTickFormatter(code="""
    return %s[tick];
    """ % (list(pos_cov_df["positions"])))
    return p

def add_annotations(bokeh_figure, transcript, max_x, max_y,
                    glyphs=["annot"], sequence=None, shift=0):
    """ Add the annotations to represent transcript annotation like boxes
    at the bottom of the plot area and sequences + codons.
    glyphs is a list and can contain: 'annot', 'nt', 'codon'.
    """
    min_y = 0
    if "annot" in glyphs:
        min_y = add_annotations_glyph(bokeh_figure, transcript, max_x, max_y)
        # annot will always be at the bottom, no test on lower_y_value
        bokeh_figure.y_range = Range1d(min_y, max_y)
    if "nt" in glyphs:
        if not sequence:
            print("Transcript Sequence is required to add nucleotides")
            raise SystemExit
        lower_y_value = add_nucleotides_glyph(bokeh_figure, max_y, sequence,
                                              shift=shift)
        if lower_y_value < min_y:
            bokeh_figure.y_range = Range1d(lower_y_value, max_y)
            min_y = lower_y_value
    if "codon" in glyphs:
        if not sequence:
            print("Transcript Sequence is required to add nucleotides")
            raise SystemExit
        lower_y_value = add_codons_glyph(bokeh_figure, max_y, sequence,
                                         shift=shift)
        if lower_y_value < min_y:
            bokeh_figure.y_range = Range1d(lower_y_value, max_y)

def add_annotations_glyph(bokeh_figure, transcript, max_x, max_y):
    """ Generate quads from annotations of the transcript's features.
    """
    dist_from_0 = max_y * 0.11
    annot_boxes_bot = dist_from_0 + (max_y * 0.05)
    margin = 0.5 # to cover the entire position at the limits

    quad_limits = {
        "left":[],
        "right":[],
        "top":[],
        "bottom":[],
        "colors":[]
    }
    for reg in ORDERED_REGIONS:
        try:
            region = getattr(transcript, reg)
        except AttributeError:
            continue
        for feature in region:
            rel_left = -1
            rel_right = -1
            if feature.start in transcript.gen_rel_map["genomic_keys"]:
                rel_left = transcript.gen_rel_map["genomic_keys"][feature.start]
            if feature.end in transcript.gen_rel_map["genomic_keys"]:
                rel_right = transcript.gen_rel_map["genomic_keys"][feature.end]

            if rel_left < 0 and rel_right < 0:
                continue # This feature is not included in the plotted area
            elif rel_left >= 0 and rel_right < 0: # cropped at the right
                rel_right = max_x + margin
            elif rel_right >= 0 and rel_left < 0: # cropped at the left
                rel_left = margin * -1
            else: # both sides of the feature are ok
                pass
            # +1 to adjust 1 based index
            rel_left -= margin
            rel_right += margin
            quad_limits["colors"].append(REG_COLORS[reg])
            # Adjut vertical
            quad_limits["top"].append(dist_from_0 * -1)
            quad_limits["bottom"].append(annot_boxes_bot * -1)
            # Adjut horizontal
            quad_limits["left"].append(rel_left) # 1 based index for annot
            quad_limits["right"].append(rel_right)

    # why on state == pre,
    source = ColumnDataSource(quad_limits)
    glyph = Quad(left="left", right="right",
                 top="top", bottom="bottom",
                 fill_color="colors")
    bokeh_figure.add_glyph(source, glyph)
    # Lower value in the Y axis, used to combined several glyphs
    lower_y_value = quad_limits["bottom"][0] - (annot_boxes_bot / 2)
    return lower_y_value

def add_nucleotides_glyph(bokeh_figure, max_y, sequence, shift=0):
    """ Generate quads for the sequence. """
    # Limits between annotation and 0 x axis line
    dist_from_0 = max_y * 0.0001
    annot_boxes_bot = dist_from_0 + (max_y * 0.02)
    box_margin = 0.5

    quad_limits = {
        "left":[],
        "right":[],
        "top":[],
        "bottom":[],
        "colors":[]
    }
    for (idx_raw, nt) in enumerate(sequence):
        idx = idx_raw + shift
        color = NT_COLORS[nt]
        quad_limits["colors"].append(color)
        quad_limits["top"].append(dist_from_0 * -1)
        quad_limits["bottom"].append(annot_boxes_bot * -1)
        quad_limits["left"].append(idx - box_margin)
        quad_limits["right"].append(idx + box_margin)

    source = ColumnDataSource(quad_limits)
    glyph = Quad(left="left", right="right",
                 top="top", bottom="bottom",
                 fill_color="colors", line_alpha=0.3,
                 line_color="white")
    lower_y_value = quad_limits["bottom"][0]
    bokeh_figure.add_glyph(source, glyph)
    return lower_y_value

def add_codons_glyph(bokeh_figure, max_y, sequence, shift=0):
    """ In a first time, only start and stop codons are displayed.
    If pre_sequence => Introns are biasing all the phasing.
    Need to remove introns from codons ?
    Or trying to introduce the phase.
    """
    dist_from_0 = max_y * 0.025
    annot_boxes_bot = dist_from_0 + (max_y * 0.025)
    margin = 1.5
    spacer = annot_boxes_bot - dist_from_0

    for frame in range(3):
        frame_color = PHASES_COLORS[frame]
        frame_seq = sequence[frame:]
        frame_codons = [frame_seq[i:i+3]
                        for i in range(0, len(frame_seq), 3)]
        # quad limits here to allow 1 line color / frame
        quad_limits = {
            "left":[],
            "right":[],
            "top":[],
            "bottom":[],
            "colors":[],
            "alpha":[],
        }

        for (idx_c, codon) in enumerate(frame_codons):
            center_nt = (idx_c * 3) + frame + 1 + shift
            try:
                color = CODONS_COLORS[codon]
                alpha = 1
            except KeyError:
                color = frame_color
                alpha = 0.1

            quad_limits["colors"].append(color)
            quad_limits["top"].append(
                (dist_from_0 * -1) - (spacer * frame))
            quad_limits["bottom"].append(
                (annot_boxes_bot * -1) - (spacer * frame))
            quad_limits["left"].append(center_nt - margin)
            quad_limits["right"].append(center_nt + margin)
            quad_limits["alpha"].append(alpha)

        source = ColumnDataSource(quad_limits)
        glyph = Quad(left="left", right="right",
                     top="top", bottom="bottom",
                     fill_color="colors", line_color=frame_color,
                     fill_alpha="alpha", line_alpha=0.3)
        lower_y_value = quad_limits["bottom"][0]
        bokeh_figure.add_glyph(source, glyph)
    return lower_y_value












####### BIC NOELYA ########
from utils import rpm_normalisation
from utils import scale_huge_peaks

def total_vbar_figure_2(pos_cov_df, title="Total VBar",
                        color="black", scale=False, rpm=False, banksize=None):
    """ Generate a figure with the total of read per position. """
    #print(pos_cov_df)
    if rpm:
        pos_cov_df[TOTAL_KEYNAME] = rpm_normalisation(list(pos_cov_df[TOTAL_KEYNAME]), banksize)
    if scale:
        pos_cov_df[TOTAL_KEYNAME] = scale_huge_peaks(list(pos_cov_df[TOTAL_KEYNAME]))
    bar_opts = dict(width=0.6, alpha=0.8)
    p = figure(title=title, plot_width=800, plot_height=350, responsive=True)
    p.vbar(bottom=0, top=pos_cov_df[TOTAL_KEYNAME].replace(0, nan),
           x=pos_cov_df.index, color=color, **bar_opts)
    p.xaxis.major_label_orientation = pi/4
    p.xaxis.formatter = FuncTickFormatter(code="""
    return %s[tick];
    """ % (list(pos_cov_df["positions"])))
    return p


def add_annotations_2(bokeh_figure, transcript, max_x, max_y,
                      glyphs=["annot"], sequence=None, shift=0, second=False):
    """ Add the annotations to represent transcript annotation like boxes
    at the bottom of the plot area and sequences + codons.
    glyphs is a list and can contain: 'annot', 'nt', 'codon'.
    """
    min_y = 0
    if "annot" in glyphs:
        min_y = add_annotations_glyph_2(bokeh_figure, transcript, max_x, max_y, second=second)
        # annot will always be at the bottom, no test on lower_y_value
        bokeh_figure.y_range = Range1d(min_y, max_y)
    if "nt" in glyphs:
        if not sequence:
            print("Transcript Sequence is required to add nucleotides")
            raise SystemExit
        lower_y_value = add_nucleotides_glyph(bokeh_figure, max_y, sequence,
                                              shift=shift)
        if lower_y_value < min_y:
            bokeh_figure.y_range = Range1d(lower_y_value, max_y)
            min_y = lower_y_value
    if "codon" in glyphs:
        if not sequence:
            print("Transcript Sequence is required to add nucleotides")
            raise SystemExit
        lower_y_value = add_codons_glyph_2(bokeh_figure, max_y, sequence,
                                         shift=shift)
        if lower_y_value < min_y:
            bokeh_figure.y_range = Range1d(lower_y_value, max_y)





def add_codons_glyph_2(bokeh_figure, max_y, sequence, shift=0):
    """ In a first time, only start and stop codons are displayed.
    If pre_sequence => Introns are biasing all the phasing.
    Need to remove introns from codons ?
    Or trying to introduce the phase.
    """
    dist_from_0 = max_y * 0.025
    annot_boxes_bot = dist_from_0 + (max_y * 0.025)
    margin = 1.5
    spacer = annot_boxes_bot - dist_from_0

    SHSH = 2
    for frame in range(3):
        frame_color = PHASES_COLORS[FRAMES_CORREL[SHSH][frame]]
        frame_seq = sequence[frame:]
        frame_codons = [frame_seq[i:i+3]
                        for i in range(0, len(frame_seq), 3)]
        # quad limits here to allow 1 line color / frame
        quad_limits = {
            "left":[],
            "right":[],
            "top":[],
            "bottom":[],
            "colors":[],
            "alpha":[],
        }

        for (idx_c, codon) in enumerate(frame_codons):
            center_nt = (idx_c * 3) + frame + 1 + shift
            try:
                color = CODONS_COLORS[codon]
                alpha = 1
            except KeyError:
                color = frame_color
                alpha = 0.1

            quad_limits["colors"].append(color)
            quad_limits["top"].append(
                (dist_from_0 * -1) - (spacer * FRAMES_CORREL[SHSH][frame]))
            quad_limits["bottom"].append(
                (annot_boxes_bot * -1) - (spacer * FRAMES_CORREL[SHSH][frame]))
            quad_limits["left"].append(center_nt - margin)
            quad_limits["right"].append(center_nt + margin)
            quad_limits["alpha"].append(alpha)

        source = ColumnDataSource(quad_limits)
        glyph = Quad(left="left", right="right",
                     top="top", bottom="bottom",
                     fill_color="colors", line_color=frame_color,
                     fill_alpha="alpha", line_alpha=0.3)
        lower_y_value = quad_limits["bottom"][0]
        bokeh_figure.add_glyph(source, glyph)
    return lower_y_value




def add_annotations_glyph_2(bokeh_figure, transcript, max_x, max_y, second=False):
    """ Generate quads from annotations of the transcript's features.
    """
    dist_from_0 = max_y * 0.11
    if second:
        dist_from_0 *= 1.5
    annot_boxes_bot = dist_from_0 + (max_y * 0.05)
    margin = 0.5 # to cover the entire position at the limits

    quad_limits = {
        "left":[],
        "right":[],
        "top":[],
        "bottom":[],
        "colors":[]
    }
    for reg in ORDERED_REGIONS:
        try:
            region = getattr(transcript, reg)
        except AttributeError:
            continue
        for feature in region:
            rel_left = -1
            rel_right = -1
            if feature.start in transcript.gen_rel_map["genomic_keys"]:
                rel_left = transcript.gen_rel_map["genomic_keys"][feature.start]
            if feature.end in transcript.gen_rel_map["genomic_keys"]:
                rel_right = transcript.gen_rel_map["genomic_keys"][feature.end]

            if rel_left < 0 and rel_right < 0:
                continue # This feature is not included in the plotted area
            elif rel_left >= 0 and rel_right < 0: # cropped at the right
                rel_right = max_x + margin
            elif rel_right >= 0 and rel_left < 0: # cropped at the left
                rel_left = margin * -1
            else: # both sides of the feature are ok
                pass
            # +1 to adjust 1 based index
            rel_left -= margin
            rel_right += margin
            quad_limits["colors"].append(REG_COLORS[reg])
            # Adjut vertical
            quad_limits["top"].append(dist_from_0 * -1)
            quad_limits["bottom"].append(annot_boxes_bot * -1)
            # Adjut horizontal
            quad_limits["left"].append(rel_left) # 1 based index for annot
            quad_limits["right"].append(rel_right)

    # why on state == pre,
    source = ColumnDataSource(quad_limits)
    glyph = Quad(left="left", right="right",
                 top="top", bottom="bottom",
                 fill_color="colors")
    bokeh_figure.add_glyph(source, glyph)
    # Lower value in the Y axis, used to combined several glyphs
    lower_y_value = quad_limits["bottom"][0] - (annot_boxes_bot / 2)
    return lower_y_value


### Ipython session ###
s = """
from loaders import Gff
gff_bic = Gff("/home/pbertin/Documents/NOELYA/LAST_NOELYA/annotations/bicistronics.gff", all_as_high=True)
gff_genes = Gff("/home/pbertin/Documents/NOELYA/LAST_NOELYA/annotations/final_at_mito.gff")
from bokeh.plotting import show, output_file
from utils import poscov_to_dataframe, get_genomic_relative_map
from bokeh_plots import total_vbar_figure_2, add_annotations, add_annotations_2
import pickle
fapick = pickle.load(open("/home/pbertin/Documents/NOELYA/LAST_NOELYA/annotations/ATMT_edite.fapick", "rb"))
import pandas as pd
bamfile = "/home/pbertin/Documents/NOELYA/LAST_NOELYA/no_pt_bams/Col0-B_norRNA_noPt_sorted.bam"
import pysam
aa = pysam.AlignmentFile(bamfile, "rb")
bs = aa.mapped
print(bs)
for f in gff_bic.all_features():
    print(f.tags)
    if "rpl5" not in f.tags["ID"]:
        continue
    tr_5p = [g.transcripts[0] for g in gff_genes.select_gene_names_id([f.tags["5p_gene"]])][0]
    tr_3p = [g.transcripts[0] for g in gff_genes.select_gene_names_id([f.tags["3p_gene"]])][0]    
    (df_5p, k1) = poscov_to_dataframe(tr_5p.mature_pos_cov(bamfile))
    (df_3p, k2) = poscov_to_dataframe(tr_3p.mature_pos_cov(bamfile))
    #print(df_5p.head())
    #print(df_3p.head())
    print(max(df_5p.positions))
    print(min(df_3p.positions))
    overlap_size = abs(max(df_5p.positions) - min(df_3p.positions))
    overlap_size = 1902
    all_df = pd.concat([df_5p, df_3p], ignore_index=True)
    all_df = all_df.groupby(['positions']).first().reset_index()
    all_df = all_df.sort_values("positions", ascending=False).reset_index()
    print(all_df.head())
    #break
    # Get relative_genomic from combination of the two transcripts    
    gen_rel = get_genomic_relative_map(all_df["positions"])
    tr_5p.gen_rel_map = gen_rel
    tr_3p.gen_rel_map = gen_rel
    fig = total_vbar_figure_2(all_df, title=f.tags["ID"], scale=True)
    add_annotations(fig, tr_5p, max(all_df.index), max(all_df["Total"]), glyphs=["annot", "codon", "nt"], sequence=tr_5p.mature_sequence(fapick))
    add_annotations_2(fig, tr_3p, max(all_df.index), max(all_df["Total"]), glyphs=["annot", "codon", "nt"], sequence=tr_3p.mature_sequence(fapick, n=(-len(tr_3p.mature_sequence(fapick)) + (overlap_size+1))), shift=len(tr_5p.mature_sequence(fapick)), second=False)
    output_file("%s_scaled.html" % (f.tags["ID"]))
    show(fig)
    #  Need to check the junction.... but seems to be closed... And then final thing, need to have a plot for a feature (to be easy to plot a region of the genome of IGORF
"""
