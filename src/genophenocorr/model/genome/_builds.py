import importlib.resources

from ._genome import Contig, GenomeBuild


def read_assembly_report(identifier: str, path: str) -> GenomeBuild:
    contigs = []
    with importlib.resources.open_text('genophenocorr.model.genome', path) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            name = fields[0]
            gb_acc = fields[4]
            refseq = fields[6]
            ucsc = fields[9]
            length = int(fields[8])
            contig = Contig(name, gb_acc, refseq, ucsc, length)
            contigs.append(contig)

    return GenomeBuild(identifier, contigs)


GRCH37 = read_assembly_report('GRCh37.p13', 'GCF_000001405.25_GRCh37.p13_assembly_report.tsv')
GRCH38 = read_assembly_report('GRCh38.p13', 'GCF_000001405.39_GRCh38.p13_assembly_report.tsv')
