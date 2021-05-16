# coding: utf-8
from loaders import Gff
from meta import MetaGene
from bokeh.plotting import show, output_file
gff_file = "/home/glihm/Documents/noelya_planchard/annotations/final_at_mito.gff"
bamfile = "/home/glihm/Documents/noelya_planchard/no_pt_bams/Col0-A_norRNA_noPt_sorted_bystrali.bam"
gff = Gff(gff_file)

trs = []
for gene in gff.genes():
    trs.append(gene.transcripts[0])
m = MetaGene(trs)
(res, fig) = m.regions_stats(bamfile, title="RPFs distribution in transcripts regions")
output_file("rpfs_distribution_regions_kmersDetails.html")
show(fig)
from loaders import Gff
from meta import MetaGene
from bokeh.plotting import show, output_file
from bokeh.layouts import gridplot
import pickle
gff_file = "/home/glihm/Documents/noelya_planchard/annotations/final_at_mito.gff"
bamfile = "/home/glihm/Documents/noelya_planchard/no_pt_bams/Col0-B_norRNA_noPt_sorted.bam"
fapick = pickle.load(open("/home/glihm/Documents/noelya_planchard/annotations/ATMT_edite.fapick", "rb"))
gff = Gff(gff_file)

stops = []
starts = []
seeked_start = ["rps4", "ccmFN2", "rps7", "nad9", "rps12", "nad4L"]
seeked_stops = ["ccmC", "nad6", "mttb", "atp1", "atp4", "cox3", "nad3"]
for gene in gff.genes():
    if gene.ID in seeked_start:
        fig = gene.get_transcript_figure(bamfile, n=200, state="cds", only_total=True, glyphs=["nt", "annot", "codon"], fapick=fapick)
        starts.append(fig)
    elif gene.ID in seeked_stops:
        fig = gene.get_transcript_figure(bamfile, n=-200, state="cds", only_total=True, glyphs=["nt", "annot", "codon"], fapick=fapick)
        stops.append(fig)
grid = gridplot([[starts[0], starts[1]], [starts[2], starts[3]], [starts[4], starts[5]]], responsive=True)
output_file("all_STARTS.html")
show(grid)
grid2 = gridplot([[stops[0], stops[1], stops[2]], [stops[3], stops[4]], [stops[5], stops[6]]], responsive=True)
output_file("all_STOPS.html")
show(grid2)
from loaders import Gff
from meta import MetaGene
from bokeh.plotting import show, output_file
import pickle
sample = "RFL23"
rep = "A"
qpcr = pickle.load(open("/home/glihm/Documents/noelya_planchard/qPCR_mean_pandas.p", "rb"))
gff_file = "/home/glihm/Documents/noelya_planchard/annotations/final_at_mito.gff"
bamfile = "/home/glihm/Documents/noelya_planchard/no_pt_bams/%s-%s_norRNA_noPt_sorted.bam" % (sample, rep)
gff = Gff(gff_file)
trs = []
for gene in gff.genes():
    trs.append(gene.transcripts[0])
m = MetaGene(trs)
(c, fig) = m.transcripts_expression(bamfile, qpcr=qpcr, sample_id=sample)
output_file("transcripts_expression_norma_qpcr_%s.html" % (sample))
show(fig)
