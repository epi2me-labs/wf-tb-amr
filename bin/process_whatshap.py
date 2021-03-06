#!/usr/bin/env python
"""Process whathap phasing info."""

import argparse
import logging

import vcf


def load_cmdline_params():
    """Return an arg parser object from arguments."""
    parser = argparse.ArgumentParser(
        description='''Process WhatsHap phasing data''')

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


def combine_phased_variants(phase_group, variants):
    """
    Combine our phased variants.

    Here we take the "lead variant", i.e. the left most and construct a new
    variant with all of those in phase ion the same codon.
    We preseve the INFO field of non-lead variants in the INFO field under
    (PREVIOUS_INFO)
    """
    refs = list()
    alts = list()
    infos = dict()
    codon = list()
    for variant in variants:

        codon.append(variant.INFO['CodonPosition'])

        refs.append(variant.REF)

        # exit if there is more than one alt - need to deal with this
        if len(variant.ALT) > 1:
            exit()

        alts.append(str(variant.ALT[0]))

        if variant.POS == phase_group:
            infos = variant.INFO
        else:
            pass

    # check for variants not in same codon
    if len(set(codon)) > 1:
        return variants

    to_remove = [
        'AA',
        'CODON_NUMBER',
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
        infos.pop(key, None)

    phased_record = vcf.model._Record(
        CHROM="NC_000962.3",
        POS=phase_group,
        ID=None,
        REF=''.join(refs),
        ALT=[vcf.model._Substitution(''.join(alts))],
        QUAL=None,
        FILTER=['PASS'],
        INFO=infos,
        sample_indexes=None,
        FORMAT=None
    )

    return [phased_record]


def process_whathap(phased_vcf, template_file, out_vcf):
    """Process whathap VCF."""
    vcf_reader = vcf.Reader(filename=phased_vcf)
    phase_groups = dict()

    template = vcf.Reader(filename=template_file)

    vcf_writer = vcf.Writer(open(out_vcf, 'w'), template)

    for record in vcf_reader:

        if is_variant_eligable(record) is False:
            vcf_writer.write_record(record)
            continue

        phase_group = check_phase_status(record)

        if phase_group is False:
            vcf_writer.write_record(record)
            continue

        if phase_group not in phase_groups:
            phase_groups[phase_group] = list()

        phase_groups[phase_group].append(record)

    # now process our phased variants
    for phase_group in phase_groups:
        phased_record = combine_phased_variants(
            phase_group, phase_groups[phase_group])
        for variant in phased_record:
            vcf_writer.write_record(variant)

    vcf_writer.close()


def main():
    """
    Squash phased variants in coding regions.

    Code to squash whatshap phased varinats in coding regions so that we
    can get the correct annotaions.
    """
    """Process our whatshap data."""
    args = load_cmdline_params().parse_args()

    logging.basicConfig(
        filename='process_whathap.log',
        filemode='w',
        format='%(asctime)s - %(message)s',
        level=logging.INFO
    )

    process_whathap(
        args.phased_vcf, args.template, args.out_vcf)


if __name__ == "__main__":
    main()
