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
    format = """##FORMAT=<ID=GT,Number=1,Type=String,
                        Description="Non-meaningful genotype">"""
    template.formats = OrderedDict(
        [vcf.parser._vcf_metadata_parser().read_format(format_string=format)]
    )
    template._column_headers.append('FORMAT')
    template._column_headers.append(sample)

    # We need to add these filters here as they use the parameters
    # set by the user at run time
    new_filters = list()

    DP_CUTOFF = args.dp
    AF_CUTOFF = args.af
    ALT_DEPTH_CUTOFF = args.alt_depth
    STRAND_BIAS_CUTOFF = args.strand_bias

    dp_header = f"""##FILTER=<ID=DP,Description="DP is less """ \
        f"""than the cut-off ({DP_CUTOFF})">"""
    af_f_header = f"""##FILTER=<ID=AF,Description="AF is less """ \
        f"""than the cut-off ({AF_CUTOFF})">"""
    ad_f_header = f"""##FILTER=<ID=ALT_DEPTH,Description="ALT Depth is """ \
        f"""less than the cut-off ({ALT_DEPTH_CUTOFF})">"""
    sb_header = f"""##FILTER=<ID=STRAND_BIAS,Description="STRAND_BIAS """ \
        f"""is greater than the cut-off """ \
        f"""({STRAND_BIAS_CUTOFF})">"""

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
    DP_CUTOFF = args.dp
    AF_CUTOFF = args.af
    ALT_DEPTH_CUTOFF = args.alt_depth
    STRAND_BIAS_CUTOFF = args.strand_bias

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

        REASON = None
        record.FILTER = list()

        ref_ADF = record.INFO['ADF'][0]
        ref_ADR = record.INFO['ADR'][0]

        for count, ALT in enumerate(record.ALT):
            new_record = copy.deepcopy(record)
            ADF = new_record.INFO['ADF'][count+1]
            ADR = new_record.INFO['ADR'][count+1]
            DP = new_record.INFO['DP']

            if DP < DP_CUTOFF:
                REASON = "DP"
                logging.info(
                    f"""\tFAILED {REASON} DP: {DP}""")
                new_record.FILTER.append('DP')
                new_record.INFO['DP'] = DP

            # calculate allele frequecny
            try:
                AF = (ADF + ADR) / DP
            except ZeroDivisionError:
                AF = 0

            # apply allele frequeny filter
            if AF < AF_CUTOFF:
                REASON = "AF"
                logging.info(
                    f"""\tFAILED {REASON} """
                    f"""ALT: {ALT} FWD: {ADF} REV: {ADR} AF: {AF}""")
                new_record.FILTER.append('AF')
                new_record.INFO['AF'] = AF

            # if total alt depth is less than a specified then apply filter
            if ADF < ALT_DEPTH_CUTOFF:
                REASON = "ALT_DEPTH"
                logging.info(
                    f"""\tFAILED {REASON} """
                    f"""ALT: {ALT} FWD: {ADF} REV: {ADR} AF: {AF}""")
                new_record.FILTER.append('ALT_DEPTH')
                new_record.INFO['ALT_DEPTH'] = ADF + ADR

            if ADR < ALT_DEPTH_CUTOFF:
                REASON = "ALT_DEPTH"
                logging.info(
                    f"""\tFAILED {REASON} """
                    f"""ALT: {ALT} FWD: {ADF} REV: {ADR} AF: {AF}""")
                new_record.FILTER.append('ALT_DEPTH')
                new_record.INFO['ALT_DEPTH'] = ADF + ADR

            # if strand bias is an issue - apply filter
            try:
                # do a fishers exact test
                d = {'forward': [ref_ADF, ADF], 'reverse': [ref_ADR, ADR]}
                table = pandas.DataFrame(data=d, index=['REF', 'ALT'])

                oddsr, p = fisher_exact(table, alternative='two-sided')
                if p == 0:
                    p = sys.float_info.min

                phred = -10 * math.log10(p)

                new_record.INFO['FS_SB'] = f'{phred:.2f}'

                if phred > STRAND_BIAS_CUTOFF:
                    logging.info(
                        f"""\tFAILED STRAND_BIAS {ALT} FWD: {ADF}/{ref_ADF}"""
                        f""" REV: {ADR}/{ref_ADR} AF: {AF}""")
                    new_record.FILTER.append('STRAND_BIAS')

            # This case we couldn't calculate strand bias
            except ZeroDivisionError:
                new_record.FILTER.append('STRAND_BIAS')
                new_record.INFO['FS_SB'] = None

            # make a VCF record - we probably need more for the info field here
            new_record.ALT = [ALT]
            new_record.INFO['AF'] = AF
            new_record.INFO['ADR'] = ADR
            new_record.INFO['ADF'] = ADF

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
