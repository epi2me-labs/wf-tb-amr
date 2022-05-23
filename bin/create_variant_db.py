#!/usr/bin/env python
"""
Take the WHO varinats excel and make a VCF.

After running:
bcftools sort `pwd`/data/primer_schemes/V2/variant_db.vcf \
    > `pwd`/data/primer_schemes/V2/variant_db.sorted.vcf
docker run -v `pwd`/data/primer_schemes/V2:/data -it zlskidmore/vt \
    vt normalize -r /data/NC_000962.3.fasta /data/variant_db.sorted.vcf \
    -o /data/variant_db.sorted.normalised.vcf
bgzip `pwd`/data/primer_schemes/V2/variant_db.sorted.normalised.vcf
tabix `pwd`/data/primer_schemes/V2/variant_db.sorted.normalised.vcf.gz
"""
import argparse
from collections import OrderedDict
import hashlib
import logging
import os
import pathlib
import urllib

from Bio import SeqIO
from Bio.Seq import Seq
from common_methods import convert_who_confidence
import numpy as np
import pandas as pd
import vcf


def load_cmdline_params():
    """Return an arg parser object from arguments."""
    who_url = "https://apps.who.int/iris/bitstream/handle/10665/341906"
    who_file = "WHO-UCN-GTB-PCI-2021.7-eng.xlsx"
    primer_schemes = "data/primer_schemes/"

    fasta = f"{primer_schemes}/V2/NC_000962.3.fasta"
    genbank = f"{primer_schemes}/V2/NC_000962.3.gb"
    template = "data/template.vcf"
    bed = f"{primer_schemes}/V2/TB_amplicons.bed"
    database = f"{primer_schemes}/V2/variant_db.vcf"

    parser = argparse.ArgumentParser(
        description='''Create VCF from WHO xlsx file''')
    parser.add_argument(
        '-r', '--reference',
        help="Reference genome in FASTA format", default=fasta)
    parser.add_argument(
        '-g', '--genbank', help="Genbank gbk file", default=genbank)
    parser.add_argument(
        '-w', '--who_variants',
        help="The WHO catalogue of MTBC mutations url.",
        default=f"{who_url}/{who_file}")
    parser.add_argument(
        '-s', '--who_sheet', help="The sheet containing the variants",
        default="Genome_indices")
    parser.add_argument(
        '-t', '--vcf_template', help="Template VCF file", default=template)
    parser.add_argument(
        '-b', '--bed_file', help="Bed file of regions to restrict to",
        default=bed)
    parser.add_argument(
        '-l', '--resistance_level',
        help="The association with resistance level to include in VCF",
        default=[1, 2, 3], nargs="+")
    parser.add_argument(
        '-e', '--group3_genes',
        help="If group 3 is included the genes to use.",
        default=['rplC', 'rrl', 'Rv0678'], nargs="+")
    parser.add_argument(
        '-i', '--git', help="Git details for VCF header")
    parser.add_argument(
        '-o', '--output_vcf_file', help="Variant database output VCF file",
        default=database)

    return parser


def translate_triplet(triplet):
    """Return a translated triplet."""
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    if triplet in table:
        return table[triplet]
    return "?"


def reverse_translate(aa):
    """Take an amino acid and get potential triplets."""
    table = {
        "I": ["ATA", "ATC", "ATT"],
        "M": ["ATG"],
        "T": ["ACA", "ACC", "ACG", "ACT"],
        "N": ["AAC", "AAT"],
        "K": ["AAA", "AAG"],
        "S": ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
        "R": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
        "L": ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
        "P": ["CCA", "CCC", "CCG", "CCT"],
        "H": ["CAC", "CAT"],
        "Q": ["CAA", "CAG"],
        "V": ["GTA", "GTC", "GTG", "GTT"],
        "A": ["GCA", "GCC", "GCG", "GCT"],
        "D": ["GAC", "GAT"],
        "E": ["GAA", "GAG"],
        "G": ["GGA", "GGC", "GGG", "GGT"],
        "F": ["TTC", "TTT"],
        "Y": ["TAC", "TAT"],
        "*": ["TAA", "TAG", "TGA"],
        "C": ["TGC", "TGT"],
        "W": ["TGG"]
    }

    if aa in table:
        return table[aa]


def get_url_data(url: str):
    """Get data from a url. We use this to retrive the WHO data."""
    user_agent = """Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7)
        Gecko/2009021910 Firefox/3.0.7"""
    headers = {'User-Agent': user_agent}
    request = urllib.request.Request(url, None, headers)
    return urllib.request.urlopen(request).read()


def get_who_data(who_variants: str, who_sheet: str) -> pd.core.frame.DataFrame:
    """Get WHO data from url or local file into a pd dataframe."""
    # Open the WHO mutations into a pandas dataframe
    if os.path.isfile(who_variants):
        data = who_variants

    if who_variants.startswith("http"):
        data = get_url_data(who_variants)

    who = pd.read_excel(data, sheet_name=who_sheet, header=[0])

    # codon number has brackets if it's an indel
    who['codon_number'] = who['codon_number'].str.strip('(')
    who['codon_number'] = who['codon_number'].str.strip(')')

    who['codon_number'] = who['codon_number'].replace(np.nan, 0)
    # get all things numeric into a numeric

    num_columns = [
        "codon_number",
        "final_annotation.Position",
        "EMB_R",
        "EMB_S",
        "INH_R",
        "INH_S",
        "PZA_R",
        "PZA_S",
        "LEV_R",
        "LEV_S",
        "RIF_R",
        "RIF_S",
        "BDQ_R",
        "BDQ_S",
        "CFZ_R",
        "CFZ_S",
        "DLM_R",
        "DLM_S",
        "AMI_R",
        "AMI_S",
        "ETH_R",
        "ETH_S",
        "KAN_R",
        "KAN_S",
        "CAP_R",
        "CAP_S",
        "MXF_R",
        "MXF_S",
        "LZD_S",
        "LZD_R"
    ]
    who[num_columns] = who[num_columns].apply(pd.to_numeric)

    who = who.fillna('.')

    return who


def process_genbank(genbank_file: str) -> list:
    """
    Process genbank.

    Process genbank to get the details we need which is a dictionary
    of genes their bases and the AA and codons.
    """
    # get all sequence records for the specified genbank file

    recs = [rec for rec in SeqIO.parse(genbank_file, "genbank")]

    # we expect only one sequence record
    if len(recs) > 1:
        raise KeyError('More than one sequence record. Expecting only one.')

    return recs[0]


def get_compliment(base: str) -> str:
    """Compliment a DNA base."""
    compliment_table = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return compliment_table[base]


def check_reference_base(position, base):
    """Check WHO reference base."""
    raise NotImplementedError


def check_reference_codon(position, base):
    """Check WHO reference codon."""
    raise NotImplementedError


def check_mutant_codon(
        genbank: list, locus_tag: str, codon_number: int, position: int,
        mutant_codon: str, alt_base: str):
    """Check a mutant codon by creating it de novo."""
    # get the refrence codon using genbank
    reference_codon = get_reference_codon_seq(
        genbank,
        locus_tag,
        codon_number
    )

    # now mutate the reference codon based on position and alt_base
    mutant_codon_check = make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )

    # make it into a format we can compare
    our_mutant_codon = ''.join([i[1] for i in mutant_codon_check])

    # do the test
    try:
        assert our_mutant_codon == mutant_codon.upper()
    except AssertionError as e:
        logging.warning(f"""CHECK_MUTANT_CODON Failed
            {locus_tag} {codon_number} {position}
            {mutant_codon} {alt_base} {mutant_codon_check} {e}""")


def check_reference_aa(features: list, codon: int, aa: str):
    """Check if the reference AA is correct."""
    # get our feature; only intrested in CDS here
    feature = [feature for feature in features if feature.type == 'CDS'][0]

    # get the amino acid from genbank
    genbank_aa = feature.qualifiers['translation'][0][codon-1]

    # do the test
    try:
        assert aa == genbank_aa
    except AssertionError as e:
        logging.warning(
            f"CHECK_REFERENCE_AA genbank says {genbank_aa} who says {aa} {e}")


def get_reference_codon_seq(
        genbank: list, locus_tag: str, codon_number: int) -> list:
    """Get a reference codon sequence based on codon number."""
    # get our feature; only intrested in CDS & rRNA here
    print(locus_tag)
    feature = [
        feature for feature in get_genbank_feature(genbank, locus_tag)
        if feature.type == 'CDS' or feature.type == 'rRNA'
    ][0]

    codon_sequence = list()

    # deal with genes on the +ve strand
    if feature.strand == 1:
        # get the genomic position nased on the codon number
        position = feature.location.start+((codon_number*3)-3)
        for i in range(position, position+3):
            codon_sequence.append((i+1, genbank.seq[i]))

    # deal with genes on the -ve strand
    elif feature.strand == -1:
        # get the genomic position nased on the codon number
        position = feature.location.end-((codon_number*3)-3)
        for i in reversed(range(position-3, position)):
            codon_sequence.append((i+1, get_compliment(genbank.seq[i])))

    return codon_sequence


def get_codon_from_ref_position(genbank, position):
    """Get a codon from a reference position."""
    pass


def make_mutant_codon(
        genbank, locus_tag: str, reference_codon: list,
        position: int, alt_base: list) -> list:
    """Given a refercence codon and a genomic mutation position and base(s)."""
    # get our feature, again onluy interest in CDS
    feature = [
        feature for feature in get_genbank_feature(genbank, locus_tag)
        if feature.type == 'CDS' or feature.type == 'rRNA'
    ][0]

    # the basis for our mutant codon is the ref codon
    mutant_codon = reference_codon.copy()

    # if +ve strand we don't need to do anything really - just make alt_base(s)
    # into a list
    if feature.strand == 1:
        alt_bases = list(alt_base)

    # if -ve strand then we need to compliment the alt_base(s)
    elif feature.strand == -1:
        alt_bases = list(Seq(alt_base).complement())

    # now for the alt_bases mutate our refrence codon codon
    for alt_base_single in alt_bases:
        mutant_codon = [
            (position, alt_base_single)
            if base[0] == position else base for base in mutant_codon
        ]
        position = position+1

    return mutant_codon


def get_genbank_feature(genbank, locus_tag):
    """Get a genbank feature given a locus tag."""
    result = list()

    # get genbank features for our locus_tag id
    for feature in genbank.features:
        if 'locus_tag' in feature.qualifiers:
            if feature.qualifiers['locus_tag'][0] == locus_tag:
                if feature.type != 'gene':
                    result.append(feature)

    return result


def prepare_vcf_header(template_vcf_file, args):
    """Add provenance to the VCF header."""
    # take our template for the header
    vcf_reader = vcf.Reader(filename=template_vcf_file)

    # add the settings used to the VCF header
    vcf_reader.metadata.copy()

    arguments = list()

    # for each of our arguments add them to a list for use in the new header
    for i in vars(args):
        if type(vars(args)[i]) is not list:
            print(type(vars(args)[i]))
            if os.path.isfile(vars(args)[i]):
                md5 = hashlib.md5(
                    pathlib.Path(vars(args)[i]).read_bytes()).hexdigest()
                arguments.append(f"<ID={i},Value={vars(args)[i]},md5={md5}>")
            else:
                arguments.append(f"<ID={i},Value={vars(args)[i]}>")

    # grab the WHO catalogue version for the file
    who_version = os.path.splitext(
            os.path.basename(vars(args)["who_variants"])
        )[0]

    # construct new metadata
    new_metadata = OrderedDict(
        origin=[
            f"<ID=WHO_MTBC_MUTATION_CATALOGUE,"
            f"Version=\"{who_version}\","
            f"Description=\"Catalogue of mutations in Mycobacterium "
            f"tuberculosis complex and their association with drug "
            f"resistance: supplementary document "
            f"https://apps.who.int/iris/handle/10665/341906\","
            f"Licence=\"CC BY-NC-SA 3.0 IGO\">"
        ],
        create_variant_db=arguments
    )

    # add our new metadat onto the old metadata
    vcf_reader.metadata.update(new_metadata)

    return vcf_reader


def load_bed_file(bed_file) -> pd.core.frame.DataFrame:
    """Load a bed file."""
    # read bed file into pandas dataframe
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    # add a meaningful header
    header = ['reference', 'start', 'end', 'gene']
    bed.set_axis(header, axis=1, inplace=True)

    return bed


def check_poistion_in_bed(
        position: int, bed: pd.core.frame.DataFrame) -> bool:
    """
    Check position.

    Check a position against a pandas dataframe of a bed file
    to see if the position is contained in any of the regions.
    """
    if len(bed[(position >= bed['start']) & (position <= bed['end'])]) > 0:
        return True
    else:
        return False


def parse_who(
        who_variants: str,
        who_sheet: str,
        template_vcf_file: str,
        output_vcf_file: str,
        genbank_file: str,
        resistance_level: int,
        bed_file: str,
        git: str,
        args):
    """Parse the WHO mutation catalogue into a useable (VCF) format."""
    # get our genbank file
    genbank = process_genbank(genbank_file)

    # prepare our vcf_template on which to base ur VCF writer on
    vcf_reader = prepare_vcf_header(template_vcf_file, args)

    # open our vcf to write to
    vcf_writer = vcf.Writer(open(output_vcf_file, 'w'), vcf_reader)

    who = get_who_data(who_variants, who_sheet)

    # check if WHO position is in the bed_file regions
    # TODO: replace this with something pandasy (the two dataframes)
    bed = load_bed_file(bed_file)
    # every row in the dataframe is a mutation in the WHO catalogue
    for index, row in who.iterrows():
        logging.info(f"PROCESSING {row['variant']}")

        pos = int(row['final_annotation.Position'])

        if check_poistion_in_bed(pos, bed) is False:
            logging.info(
                f"""DISCARDING: {pos} is NOT within our defined regions.""")
            continue
        else:
            logging.info(
                f"""RETAINING: {pos} is within our defined regions.""")

        # only interested in those that have a resistance confidence rating
        # as input by user
        conf_grades = [
            int(convert_who_confidence(value)) for key, value in
            row.iteritems() if key.endswith("Conf_Grade")]

        # if resistance_level is in confidence grades
        discard = None
        for res in resistance_level:
            if int(res) in conf_grades:
                if int(res) != 3:
                    logging.info(
                        f"""RETAINING: Requested {res} \
                        is in {conf_grades}.""")
                    discard = False
                else:
                    if row['final_annotation.Gene'] in args.group3_genes:
                        logging.info(
                            f"""RETAINING: Requested {res} \
                            is in {conf_grades}.""")
                        discard = False
            else:
                logging.info(
                    f"""DISCARDING: Requested {res} \
                    is NOT in {conf_grades}.""")
                if discard is not False:
                    discard = True

        if discard is True:
            continue

        logging.info("VARIANT RETAINED")

        # collect all ratings for all antibiotics for the VCF
        antibiotic_confidence = [
            f"{key.replace('_Conf_Grade', '')}|{convert_who_confidence(value)}"
            for key, value in row.iteritems() if key.endswith("Conf_Grade")]

        # get the genbank features for our gene
        features = get_genbank_feature(
            genbank, row['final_annotation.LocusTag'])

        # in the WHO database anything with a codon number
        # is in a coding region
        if row['codon_number'] != 0:

            reference_codon_seq = get_reference_codon_seq(
                genbank, row['final_annotation.LocusTag'], row['codon_number'])

            check_mutant_codon(
                genbank,
                row['final_annotation.LocusTag'],
                row['codon_number'],
                row['final_annotation.Position'],
                row['alt_nt'],
                row['final_annotation.AlternativeNucleotide'].upper()
            )

            make_mutant_codon(
                genbank,
                row['final_annotation.LocusTag'],
                reference_codon_seq,
                row['final_annotation.Position'],
                row['final_annotation.AlternativeNucleotide'].upper())

        # make our vcf record
        record = vcf.model._Record(
            CHROM=row['final_annotation.Reference'],
            POS=row['final_annotation.Position'],
            ID=row['variant'],
            REF=row['final_annotation.ReferenceNucleotide'].upper(),
            ALT=[vcf.model._Substitution(
                    row['final_annotation.AlternativeNucleotide'].upper())],
            QUAL=None,
            FILTER=None,
            INFO=dict(
                # INDEL=True,
                CODON_NUMBER=row['codon_number'],
                AA=None,
                EFFECT=row['final_annotation.Effect'],
                GENE=row['final_annotation.Gene'],
                GENE_LOCUS=row['final_annotation.LocusTag'],
                HGVS_NUCLEOTIDE=row[
                    'final_annotation.TentativeHGVSNucleotidicAnnotation'],
                HGVS_PROTEIN=row[
                    'final_annotation.TentativeHGVSProteicAnnotation'],
                PROTEIN_ID=row['final_annotation.ProteinId'],
                STRAND=[
                    feature.strand for feature in get_genbank_feature(
                        genbank, row['final_annotation.LocusTag'])
                    if feature.type != 'gene'
                ],
                ANTIBIOTICS=','.join(antibiotic_confidence),
                FEATURE_TYPE=','.join([feature.type for feature in features]),
                ORIGIN='WHO_CANONICAL'
            ),
            FORMAT=None,
            sample_indexes=None
        )
        print(record)
        vcf_writer.write_record(record)
    vcf_writer.close()


def main():
    """Create a VCF variant database from the WHO variants."""
    logging.basicConfig(
        # filename='create_variant_db.log',
        # filemode='w',
        format='%(asctime)s - %(message)s',
        level=logging.INFO
    )

    args = load_cmdline_params().parse_args()

    parse_who(
        args.who_variants,
        args.who_sheet,
        args.vcf_template,
        args.output_vcf_file,
        args.genbank,
        args.resistance_level,
        args.bed_file,
        args.git,
        args
    )


if __name__ == "__main__":
    main()
