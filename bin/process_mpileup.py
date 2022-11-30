#!/usr/bin/env python
"""Return variants from a bam file."""

import argparse
from collections import namedtuple, OrderedDict
import copy
import logging
import math
import sys

import pandas
from scipy.stats import fisher_exact
import vcf


def load_cmdline_params():
    """Return an arg parser object from arguments."""
    parser = argparse.ArgumentParser(
        description='''Process BCFTools mpileup data''')
    parser.add_argument(
        '-t', '--template', required=True, dest="template",
        help="VCF tenplate file")
    parser.add_argument(
        '-m', '--mpileup', required=True, dest="mpileup",
        help="BCFTools mpileup VCF file")
    parser.add_argument(
        '-s', '--sample', required=True, dest="sample",
        help="Sample name")
    parser.add_argument(
        '-a', '--af', required=True, type=float, default=0.1,
        help="MAF cutoff", dest="af")
    parser.add_argument(
        '-p', '--dp', required=True, type=float, default=20,
        help="DP cutoff", dest="dp")
    parser.add_argument(
        '-d', '--alt_depth', required=True, type=int, default=10,
        help="ALT depth cutoff", dest="alt_depth")
    parser.add_argument(
        '-b', '--strand_bias', required=True, type=float, default=150,
        help="Strand bias phred cutoff", dest="strand_bias")
    parser.add_argument(
        '-o', '--out_vcf', required=True, dest="out_vcf",
        help="Processed VCF file")

    return parser


def setup_template(template, sample, args):
    """Add required VCF headers."""
    template = vcf.Reader(filename=template)

    # add required stuff to the universal template
    var_format = """##FORMAT=<ID=GT,Number=1,Type=String,
                        Description="Non-meaningful genotype">"""
    template.formats = OrderedDict(
        [vcf.parser._vcf_metadata_parser().read_format(
            format_string=var_format)]
    )
    template._column_headers.append('FORMAT')
    template._column_headers.append(sample)

    # We need to add these filters here as they use the parameters
    # set by the user at run time
    new_filters = list()

    dp_cutoff = args.dp
    ad_cutoff = args.af
    alt_depth_cutoff = args.alt_depth
    strand_bias_cutoff = args.strand_bias

    dp_header = f"""##FILTER=<ID=DP,Description="DP is less """ \
        f"""than the cut-off ({dp_cutoff})">"""
    af_f_header = f"""##FILTER=<ID=AF,Description="AF is less """ \
        f"""than the cut-off ({ad_cutoff})">"""
    ad_f_header = f"""##FILTER=<ID=ALT_DEPTH,Description="ALT Depth is """ \
        f"""less than the cut-off ({alt_depth_cutoff})">"""
    sb_header = f"""##FILTER=<ID=STRAND_BIAS,Description="STRAND_BIAS """ \
        f"""is greater than the cut-off """ \
        f"""({strand_bias_cutoff})">"""

    new_filters.append(
        vcf.parser._vcf_metadata_parser().read_filter(dp_header))
    new_filters.append(
        vcf.parser._vcf_metadata_parser().read_filter(af_f_header))
    new_filters.append(
        vcf.parser._vcf_metadata_parser().read_filter(ad_f_header))
    new_filters.append(
        vcf.parser._vcf_metadata_parser().read_filter(sb_header))

    template.filters = OrderedDict(new_filters)

    return template


def process_mpileup(mpileup, template, out_vcf, sample, args):
    """Process mpileup."""
    vcf_reader = vcf.Reader(filename=mpileup)

    # open our universal template VCF file
    template = setup_template(template, sample, args)

    # open our output VCF for writing and base it on our universal template
    vcf_writer = vcf.Writer(open(out_vcf, 'w'), template)
    processed_records = list()

    # set our cutoffs
    dp_cutoff = args.dp
    ad_cutoff = args.af
    alt_depth_cutoff = args.alt_depth
    strand_bias_cutoff = args.strand_bias

    # for every record in our pileup decide if a variant and then write
    # to whatshap vcf for phasing
    for record in vcf_reader:

        logging.info(f"""CONSIDERING: {record.POS}
                         REF: {record.REF}
                         FWD: {record.INFO['ADF'][0]}
                         REV: {record.INFO['ADR'][0]}""")

        # those that have not been annotated are not in our database
        # so - we need to keep these in case of MNPs
        # if 'FEATURE_TYPE' not in record.INFO:
        #     logging.info("\tNOT IN DATABASE")
        #     continue

        reason = None
        record.FILTER = list()

        ref_adf = record.INFO['ADF'][0]
        ref_adr = record.INFO['ADR'][0]

        for count, alt in enumerate(record.ALT):
            new_record = copy.deepcopy(record)
            adf = new_record.INFO['ADF'][count+1]
            adr = new_record.INFO['ADR'][count+1]
            dp = new_record.INFO['DP']

            if dp < dp_cutoff:
                reason = "DP"
                logging.info(
                    f"""\tFAILED {reason} DP: {dp}""")
                new_record.FILTER.append('DP')
                new_record.INFO['DP'] = dp

            # calculate allele frequecny
            try:
                af = (adf + adr) / dp
            except ZeroDivisionError:
                af = 0

            # apply allele frequeny filter
            if af < ad_cutoff:
                reason = "AF"
                logging.info(
                    f"""\tFAILED {reason} """
                    f"""ALT: {alt} FWD: {adf} REV: {adr} AF: {af}""")
                new_record.FILTER.append('AF')
                new_record.INFO['AF'] = af

            # if total alt depth is less than a specified then apply filter
            if adf < alt_depth_cutoff:
                reason = "ALT_DEPTH"
                logging.info(
                    f"""\tFAILED {reason} """
                    f"""ALT: {alt} FWD: {adf} REV: {adr} AF: {af}""")
                new_record.FILTER.append('ALT_DEPTH')
                new_record.INFO['ALT_DEPTH'] = adf + adr

            if adr < alt_depth_cutoff:
                reason = "ALT_DEPTH"
                logging.info(
                    f"""\tFAILED {reason} """
                    f"""ALT: {alt} FWD: {adf} REV: {adr} AF: {af}""")
                new_record.FILTER.append('ALT_DEPTH')
                new_record.INFO['ALT_DEPTH'] = adf + adr

            # if strand bias is an issue - apply filter
            if 'IGNORE_SB' not in record.INFO:

                try:
                    # do a fishers exact test
                    d = {'forward': [ref_adf, adf], 'reverse': [ref_adr, adr]}
                    table = pandas.DataFrame(data=d, index=['REF', 'ALT'])

                    oddsr, p = fisher_exact(table, alternative='two-sided')
                    if p == 0:
                        p = sys.float_info.min

                    phred = -10 * math.log10(p)

                    new_record.INFO['FS_SB'] = f'{phred:.2f}'

                    if phred > strand_bias_cutoff:
                        logging.info(
                            f"\tFAILED STRAND_BIAS {alt} FWD: {adf}/{ref_adf}"
                            f" REV: {adr}/{ref_adr} AF: {af}")
                        new_record.FILTER.append('STRAND_BIAS')

                # This case we couldn't calculate strand bias
                except ZeroDivisionError:
                    new_record.FILTER.append('STRAND_BIAS')
                    new_record.INFO['FS_SB'] = None

            else:
                new_record.INFO['FS_SB'] = None

            # make a VCF record - we probably need more for the info field here
            new_record.ALT = [alt]
            new_record.INFO['AF'] = af
            new_record.INFO['ADR'] = adr
            new_record.INFO['ADF'] = adf

            processed_records.append(new_record)

    for record in processed_records:

        new_record = vcf.model._Record(
                CHROM=record.CHROM,
                POS=record.POS,
                ID=None,
                REF=record.REF,
                ALT=record.ALT,
                QUAL=None,
                FILTER=record.FILTER,
                INFO=record.INFO,
                FORMAT='GT',
                sample_indexes=[sample]
            )
        calldata = namedtuple('CallData', ['GT'])
        sample = vcf.model._Call(new_record, sample, data=calldata(GT='0/1'))
        new_record.samples = [sample]

        vcf_writer.write_record(new_record)

    vcf_writer.close()


def main():
    """Process mpileup data."""
    logging.basicConfig(
        filename='process_mpileup.log',
        filemode='w',
        format='%(asctime)s - %(message)s',
        level=logging.INFO
    )

    args = load_cmdline_params().parse_args()

    process_mpileup(
        args.mpileup, args.template, args.out_vcf, args.sample, args
    )


if __name__ == "__main__":
    main()
