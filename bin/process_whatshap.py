#!/usr/bin/env python
"""Process whathap phasing info."""

import argparse
from collections import namedtuple, OrderedDict
import logging

from Bio.Seq import Seq
import vcf


def load_cmdline_params():
    """Return an arg parser object from arguments."""
    parser = argparse.ArgumentParser(
        description='''Process WhatsHap phasing data''')

    parser.add_argument(
        '-s', '--sample', required=True, dest="sample",
        help="Sample identifier")
    parser.add_argument(
        '-p', '--phased_vcf', required=True, dest="phased_vcf",
        help="whatshap phased VCF")
    parser.add_argument(
        '-t', '--template', required=True, dest="template",
        help="VCF template file")
    parser.add_argument(
        '-o', '--out_vcf', required=True, dest="out_vcf",
        help="Processed VCF file")

    return parser


def make_call(record, sample):
    """Give our VCF a sample to make it compatible with downstream stuff."""
    data = namedtuple("calldata", ["GT"])
    sample = vcf.model._Call(record, sample, data('1/1'))
    record.FORMAT = 'GT'
    record.samples = [sample]
    return record


def is_variant_eligable(record):
    """Check if a variant is in CDS coding region."""
    if record.INFO['CodonPosition'] is None:
        return False

    if int(record.INFO['CodonPosition']) > 0:
        return True

    return False


def check_phase_status(record):
    """Check if a variant is phased."""
    for call in record.samples:

        if 'PS' in record.FORMAT:
            phase_group = getattr(call.data, 'PS')
            logging.info(f"""
                PHASE GROUP: {phase_group}""")
            return phase_group

    return False


def check_codon_status(sample, phase_group, vcf_writer):
    """
    Check the codon our phased variants is in.

    We need to know if they are in the same codon or not. If not just write to
    file. If yes keep for processing.
    """
    new_phase_groups = dict()
    for position in phase_group:
        codons = dict()
        for variant in phase_group[position]:

            if variant.INFO['CodonPosition'] not in codons:
                codons[variant.INFO['CodonPosition']] = list()
            codons[variant.INFO['CodonPosition']].append(variant)

        for codon in codons:
            if len(codons[codon]) > 1:
                if position not in new_phase_groups:
                    new_phase_groups[position] = dict()
                new_phase_groups[position][codon] = codons[codon]
            else:
                record = make_call(codons[codon][0], sample)
                vcf_writer.write_record(record)

    return new_phase_groups


def unique(sequence):
    """Get uniques in lists for info."""
    # sometimes info items are lists, join if so
    new_sequence = list()

    for part in sequence:
        if type(part) == list:
            part = "|".join(filter(None, part))
        new_sequence.append(part)

    # retrun only unique values
    seen = set()
    return [x for x in new_sequence if not (x in seen or seen.add(x))]


def combine_phased_variants(sample, phase_group, codon, variants):
    """
    Combine our phased variants.

    Here we take the "lead variant", i.e. the left most and construct a new
    variant with all of those in phased in the same codon.
    We preseve the INFO field of non-lead variants in the INFO field under
    PREVIOUS_INFO
    """
    refs = list()
    alts = list()

    codon_positions = {0: None, 1: None, 2: None}

    # sort the variants in the phase group in position order
    sorted_variants = sorted(variants, key=lambda x: x.POS, reverse=False)
    info_items = sorted_variants[0].INFO
    ref_codon = sorted_variants[0].INFO['RefCodon']
    if type(ref_codon) == list:
        ref_codon = ref_codon[0]
    # if on negative strand need to reverse complement codon
    codon_lookup = {0: 2, 1: 1, 2: 0}

    try:
        strand, mut = sorted_variants[0].INFO['Comments'][0].split(":")
    except AttributeError:
        strand = "Positive"

    if strand == "Negative":
        dna = Seq(ref_codon)
        ref_codon = dna.reverse_complement()

    # assign variants in phase group to a position in the codon
    for variant in sorted_variants:
        codon_position = variant.INFO['SNPCodonPosition']
        if strand == "Negative":
            # if on negative strand need to reverse codon positions
            codon_position = codon_lookup[codon_position]
        codon_positions[codon_position] = variant

    new_infos = dict()

    # process each codon position
    for codon in codon_positions:
        variant = codon_positions[codon]

        # a variant might not be present at each codon position
        if variant is not None:

            # we don't handle multiple alts so exit
            if len(variant.ALT) > 1:
                exit()

            # append the ref and the alt for the variant
            refs.append(variant.REF)
            alts.append(str(variant.ALT[0]))

            # add all info items so we can use them for the MNP
            for key in info_items:
                if key not in new_infos:
                    new_infos[key] = list()
                new_infos[key].append(variant.INFO[key])

        # if no variant present in the middle of the codon then we need to pad
        # the reference and the alt
        if variant is None and codon == 1:

            refs.append(ref_codon[codon])
            alts.append(ref_codon[codon])

            # add None to the infos for this position
            for key in info_items:
                if key not in new_infos:
                    new_infos[key] = list()
                new_infos[key].append(None)

    # combine all infos that have duplicate items
    to_remove = [
        'AA',
        'CODON_NUMBER',
        'I16',
        'QS',
        'AD',
        'ADF',
        'ADR',
        'ANTIBIOTICS',
        'GENE',
        'GENE_LOCUS',
        'HGVS_NUCLEOTIDE',
        'HGVS_PROTEIN',
        'ORIGIN',
        'FEATURE_TYPE',
        'STRAND',
        'EFFECT',
        'PROTEIN_ID'
    ]

    for key in to_remove:
        new_infos.pop(key, None)

    # we have to now combine info items
    final_infos = dict()

    for key in new_infos:
        unique_items = unique(new_infos[key])

        if len(unique_items) > 1:
            final_infos[key] = unique_items
        else:
            final_infos[key] = unique_items[0]

    # this is our new VCF record
    phased_record = vcf.model._Record(
        CHROM="NC_000962.3",
        POS=sorted_variants[0].POS,
        ID=None,
        REF=''.join(refs),
        ALT=[vcf.model._Substitution(''.join(alts))],
        QUAL=None,
        FILTER=['PASS'],
        INFO=final_infos,
        sample_indexes=None,
        FORMAT='GT'
    )

    phased_record = make_call(phased_record, sample)
    return [phased_record]


def process_whathap(sample, phased_vcf, template_file, out_vcf):
    """Process whathap VCF."""
    vcf_reader = vcf.Reader(filename=phased_vcf)
    phase_groups = dict()

    template = vcf.Reader(filename=template_file)

    # add required stuff to the universal template
    format_field = """##FORMAT=<ID=GT,Number=1,Type=String,
                        Description="Non-meaningful genotype">"""
    template.formats = OrderedDict(
        [vcf.parser._vcf_metadata_parser().read_format(
            format_string=format_field)]
    )
    template._column_headers.append('FORMAT')
    template._column_headers.append(sample)

    vcf_writer = vcf.Writer(open(out_vcf, 'w'), template)

    for record in vcf_reader:

        # check if a variant is coding or not - won't affect annottaion
        if is_variant_eligable(record) is False:
            record = make_call(record, sample)
            vcf_writer.write_record(record)
            continue

        # check if a variant is phased
        phase_group = check_phase_status(record)

        if phase_group is False:
            record = make_call(record, sample)
            vcf_writer.write_record(record)
            continue

        # add variant to phase groups
        if phase_group not in phase_groups:
            phase_groups[phase_group] = list()

        phase_groups[phase_group].append(record)

    # for our phased varinats check if they are in the same codon and make
    # new groups if they are
    codon_aware_phase_groups = check_codon_status(
        sample, phase_groups, vcf_writer)

    # now process our codon aware phased variants
    for phase_group in codon_aware_phase_groups:
        for codon in codon_aware_phase_groups[phase_group]:
            phased_record = combine_phased_variants(
                sample,
                phase_group,
                codon,
                codon_aware_phase_groups[phase_group][codon])
            for variant in phased_record:
                record = make_call(record, sample)
                vcf_writer.write_record(variant)

    vcf_writer.close()


def main():
    """
    Squash phased variants in coding regions.

    Code to squash whatshap phased variants in coding regions so that we
    can get the correct annotaions.
    """
    args = load_cmdline_params().parse_args()

    logging.basicConfig(
        filename='process_whathap.log',
        filemode='w',
        format='%(asctime)s - %(message)s',
        level=logging.INFO
    )

    process_whathap(
        args.sample, args.phased_vcf, args.template, args.out_vcf)


if __name__ == "__main__":
    main()
