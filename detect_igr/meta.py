#!/usr/bin/python3

""" Meta module handles the data from a pool of genes. """

# Stdlib
from collections import defaultdict
import os

# Third
from numpy import mean, median

# Local
from constants import *
from sambam import get_stranded_region_coverage_reduced
from sambam import get_bank_size
from utils import interval_around
#from pygal_plots import regions_stats_plot
from bokeh_plots import regions_stats_plot
from bokeh_plots import frames_plot
from bokeh_plots import envs_plot
from bokeh_plots import plot_expression

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"


class MetaGene():
    """ MetaGenes base class. """
    def __init__(self, transcripts):
        """ Init with a list of transcripts. """
        self.transcripts = transcripts

    # Private methods
    def __repr__(self):
        return "<%s.%s, #transcripts=%s>" % (__name__, self.__class__.__name__,
                                             len(self.transcripts))

    # Public methods
    def regions_stats(self, bamfile, shift=SHIFT, side=SIDE, title="Regions",
                      only_total=False, colors=["black", "gray", "gray"]):
        """ One loop over all transcripts to have the
        total number of reads on each region + the sizes.
        """
        regs_stats = {"counts":{}, "sizes":{}}
        for reg_name in ORDERED_REGIONS:
            regs_stats["counts"][reg_name] = {}
            regs_stats["sizes"][reg_name] = []
        for transcript in self.transcripts:
            # Get all counts for regions + the size
            for reg_name in ORDERED_REGIONS:
                try:
                    region = getattr(transcript, reg_name)
                    reg_size = region.length()
                    reg_counts = region.counts(bamfile, shift=shift, side=side)
                except AttributeError:
                    continue
                regs_stats["sizes"][reg_name].append(reg_size)
                for kmer in reg_counts:
                    try:
                        regs_stats["counts"][reg_name][kmer] += reg_counts[kmer]
                    except KeyError:
                        regs_stats["counts"][reg_name][kmer] = reg_counts[kmer]
        bank_size = get_bank_size(bamfile)
        # Each kmer count is divided by the mean of regions sizes
        for reg in regs_stats["counts"]:
            reg_mean = mean(regs_stats["sizes"][reg])
            for (kmer, cov) in regs_stats["counts"][reg].items():
                regs_stats["counts"][reg][kmer] = round(cov / reg_mean, 1)

        regs_figure = regions_stats_plot(regs_stats, bamfile,
                                         only_total=only_total,
                                         colors=colors,
                                         title=title)
        return (regs_stats, regs_figure)

    def allframes_counts(self, bamfile, shift=SHIFT, side=SIDE,
                         only_total=False, title="All frames"):
        """ Return the pooled counts in each frame of the CDSes
        of the transcripts.
        """
        frames_counts = {0:{}, 1:{}, 2:{}}

        for transcript in self.transcripts:
            counts = transcript.CDS.allframes_reads(
                bamfile, shift=shift, side=side)
            for frame in counts:
                for kmer in counts[frame]:
                    kmer_count = counts[frame][kmer]
                    try:
                        frames_counts[frame][kmer] += kmer_count
                    except KeyError:
                        frames_counts[frame][kmer] = kmer_count
        bamname = os.path.basename(bamfile).split(".")[0]
        title = "frames counts {bamname} (shift={shift}, side={side}')".format(
            **locals())
        frames_fig = frames_plot(frames_counts, bamfile,
                                 only_total=only_total, title=title)
        return (frames_counts, frames_fig)

    def start_stop_envs(self, bamfile, shift=SHIFT, side=SIDE, kmers=KMERS,
                        inner_interval=INNER, outer_interval=OUTER,
                        only_total=False):
        """ Return a dictionary with the start and stop meta coverage. """
        meta_start = {}
        meta_stop = {}
        default_required = True # Need to catch all kmers => default dic int

        if kmers != [TOTAL_KEYNAME]: # Add the total key if kmers selected
            kmers.append(TOTAL_KEYNAME)
            default_required = False

        for i, j in zip(range(outer_interval * -1, inner_interval + 1),
                        range(inner_interval * -1, outer_interval + 1)):
            # Enforce the selection of the selected kmers
            if not default_required:
                meta_start[i] = {}.fromkeys(kmers, 0)
                meta_stop[j] = {}.fromkeys(kmers, 0)
            else:
                meta_start[i] = defaultdict(int)
                meta_start[i][TOTAL_KEYNAME] = 0
                meta_stop[j] = defaultdict(int)
                meta_stop[j][TOTAL_KEYNAME] = 0

        for transcript in self.transcripts:
            transcript.CDS.start_around(bamfile,
                                        meta_start,
                                        shift=shift,
                                        side=side,
                                        inner=inner_interval,
                                        outer=outer_interval,
                                        kmers=kmers)
            transcript.CDS.stop_around(bamfile,
                                       meta_stop,
                                       shift=shift,
                                       side=side,
                                       inner=inner_interval,
                                       outer=outer_interval,
                                       kmers=kmers)
        envs_fig = envs_plot(meta_start, meta_stop, bamfile,
                             only_total=only_total)
        return (meta_start, meta_stop, envs_fig)

    def transcripts_expression(self, bamfile, shift=SHIFT, side=SIDE,
                               kmers=KMERS, sample_id=None,
                               rnaseq=None, qpcr=None,
                               title="Expression profile"):
        """ compute the expression profile of the transcripts, normalising
        by rnaseq counts or qpcr abundances if provided. rnaseq / qpcr are
        exclusive. They must be pandas dataframes with genes id / name as
        index.
        """
        counts = {}
        for tr in self.transcripts:
            tr_counts = tr.counts(
                bamfile, shift=shift, side=side, kmers=kmers)
            tr_total = sum(
                [tr_counts[k][TOTAL_KEYNAME]
                 for k in tr_counts if tr_counts[k]])
            if not qpcr.empty:
                gene_id = tr.get_tag("gene_id")
                tr_total = tr_total / qpcr.loc[gene_id][sample_id]
            counts[gene_id] = tr_total

            types = ["Complex V",
                     "Complex V",
                     "Complex V",
                     "Complex V",
                     "Complex V",
                     "Complex V",
                     "C-type cytochrome maturation",
                     "C-type cytochrome maturation",
                     "C-type cytochrome maturation",
                     "C-type cytochrome maturation",
                     "C-type cytochrome maturation",
                     "Complex III",
                     "Complex IV",
                     "Complex IV",
                     "Complex IV",
                     "Maturase",
                     "Protein Transport",
                     "Complex I",
                     "Complex I",
                     "Complex I",
                     "Complex I",
                     "Complex I",
                     "Complex I",
                     "Complex I",
                     "Complex I",
                     "Complex I",
                     "Ribosomal Proteins",
                     "Ribosomal Proteins",
                     "Ribosomal Proteins",
                     "Ribosomal Proteins",
                     "Ribosomal Proteins",
                     "Ribosomal Proteins",
                     "Ribosomal Proteins"]

        if not title:
            title_ok = "Expression profile: %s" % (sample_id)
        else:
            title_ok = title
        fig = plot_expression(counts, types=types, title=title_ok)

        return (counts, fig)
