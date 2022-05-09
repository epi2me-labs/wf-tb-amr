#!/usr/bin/env python
"""Methods common to reporting."""

import pandas as pd
import vcf


def call_resistance(resistance: dict, antibiotics: dict) -> str:
    """Call the antibiotic resistance from variant data."""
    if len(resistance['resistant'].keys()) == 0:
        resistance['resistance_level'] = "NONE"

    if len(resistance['resistant'].keys()) > 0:
        resistance['resistance_level'] = 'RES'

    fluoroquinolones = [
        antibiotic for antibiotic in antibiotics
        if antibiotics[antibiotic]['drug-class'] == 'fluoroquinolones']

    groupa = [
        antibiotic for antibiotic in antibiotics
        if antibiotics[antibiotic]['group'] == 'A']

    res = resistance['resistant']
    # MDR-TB: Resistance of Mycobacterium tuberculosis strains to at
    # least isoniazid and rifampicin, the cornerstone medicines for the
    # treatment of TB. Rifampicin-resistant disease on its own requires
    # similar clinical management as MDR-TB.

    if 'RIF' in res:
        resistance['resistance_level'] = 'RR'

        if 'RIF' in res and 'INH' in res:
            resistance['resistance_level'] = 'MDR'

    if resistance['resistance_level'] in ['MDR', 'RR']:

        fluro_count = 0

        groupa_count = 0

        fluro_groupa_count = 0

        for fluro in fluoroquinolones:
            if fluro in res:
                fluro_count += 1
                # remove those fluros we've already found from group a
                if fluro in groupa:
                    groupa.remove(fluro)
                    fluro_groupa_count += 1

        for a in groupa:
            if a in res:
                groupa_count += 1

        # The new definition of pre-XDR-TB is: TB caused by Mycobacterium
        # tuberculosis (M. tuberculosis)strains that fulfil the definition of
        # multidrug resistant and rifampicin-resistant TB (MDR/RR-TB) and
        # which are also resistant to any fluoroquinolone.

        if fluro_count >= 1:
            resistance['resistance_level'] = 'pre-XDR'

        # The updated definition of XDR-TB is: TB caused by Mycobacterium
        # tuberculosis (M. tuberculosis)strains that fulfil the definition
        # of MDR/RR-TB and which are also resistant to any fluoroquinolone
        # and at least one additional Group A drug (Group A drugs are the
        # most potent group of drugs in the ranking of second-line medicines
        # for the treatment of drug-resistant forms of TB using
        # longer treatment regimens and comprise levofloxacin, moxifloxacin,
        # bedaquiline and linezolid).

        # This is the case where it's a fluro and a non fluro group a
        if fluro_count >= 1 and groupa_count >= 1:
            resistance['resistance_level'] = 'XDR'

        # This is the case where it's a flouro and a fluro group a
        if fluro_count >= 2 and fluro_groupa_count >= 1:
            resistance['resistance_level'] = 'XDR'

    return resistance


def process_resistance(
        vcf_file: dict, antibiotics: dict, set_score: int) -> dict:
    """Calculate resistance level - None, MDR or XDR."""
    result = dict(
        resistance_level=None,
        resistant=dict(),
        susceptible=dict()
    )

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        for antibiotic_score in record.INFO['ANTIBIOTICS']:
            antibiotic, score = antibiotic_score.split('|')

            if int(score) == set_score:
                if antibiotic not in result['resistant']:
                    result['resistant'][antibiotic] = \
                        antibiotics[antibiotic]
                    result['resistant'][antibiotic]['variants'] = list()
                result['resistant'][antibiotic]['variants'].append(record)

    for antibiotic in antibiotics:
        if antibiotic not in result['resistant']:
            if antibiotics[antibiotic]['in-assay'] == "True":
                result['susceptible'][antibiotic] = antibiotics[antibiotic]
                result['susceptible'][antibiotic]['variants'] = list()

    return result


def comma_separator(sequence: list) -> str:
    """Grammatically correct comma separated list printer."""
    if not sequence:
        return ''
    if len(sequence) == 1:
        return sequence[0]
    return '{} and {}'.format(', '.join(sequence[:-1]), sequence[-1])


def thresh_func(coverage_thresh, x):
    """Determine if coverage meets a pre-set threshold."""
    if int(coverage_thresh) < 0:
        return (x >= -int(coverage_thresh)).sum()/len(x)

    return (x >= int(coverage_thresh)).sum()/len(x)


def process_coverage(
        coverage_file: str, threshold: str) -> dict:
    """Process coverage data for a sample."""
    coverage = pd.read_csv(coverage_file, sep='\t', header=None)

    column_names = ['ref', 'start', 'end', 'gene', 'position', 'coverage']
    coverage = coverage.set_axis(column_names, axis=1, inplace=False)

    # add thresholds to our coverage_groupings
    coverage_thresh, target_thresh = threshold.split(',')

    summary = coverage.groupby('gene')['coverage'].agg(
        mean=pd.DataFrame.mean,
        median=pd.DataFrame.median,
        threshold=lambda x: thresh_func(coverage_thresh, x),
        lt_20x_bases=lambda x: (x < 20).sum(),
        ge_20x_bases=lambda x: (x >= 20).sum(),
        ge_100x_bases=lambda x: (x >= 100).sum(),
        lt_20x_percent=lambda x: (x < 20).sum()/len(x),
        ge_20x_percent=lambda x: (x >= 20).sum()/len(x),
        ge_100x_percent=lambda x: (x >= 100).sum()/len(x))

    return summary


def determine_status(coverage_summary: pd.DataFrame, threshold: str) -> dict:
    """
    Determine amplicon coverage status.

    positive: -20,2 - fail is less than 20 reads in more than 2 amplicons
    ntc: 20,3 - fail is more than 20 reads in more than 3 amplicons
    sample: -20,11 - fail is less than 20 reads in more than 11 amplicons
    """
    result = dict()

    coverage, targets = threshold.split(',')

    coverage = int(coverage)
    targets = int(targets)

    # coverage is negative which means fail is less than that number
    if coverage < 0:
        coverage_summary['passed'] = (coverage_summary['median'] > -coverage)

    # coverage is positive which means fail is more than that number
    if coverage > 0:
        coverage_summary['passed'] = (coverage_summary['median'] < coverage)

    failed_count = coverage_summary[~coverage_summary['passed']].shape[0]

    result['status'] = 'pass'

    # targets is positive which means fail is more than that number
    if targets > 0:
        if failed_count > targets:
            result['status'] = 'fail'

    # targets is negative which means fail is less than that number
    if targets < 0:
        if failed_count < -targets:
            result['status'] = 'fail'

    result['failed_targets'] = coverage_summary[
        ~coverage_summary['passed']].to_dict('index')

    result['passed_targets'] = coverage_summary[
        coverage_summary['passed']].to_dict('index')

    return result


def variants_table_from_vcf(
        vcf_file: str, info_fields: list, confidence: int) -> pd.DataFrame:
    """Give a VCF file return a pandas data frame."""
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    data = dict(
        chromosome=list(),
        position=list(),
        reference=list(),
        alternate=list()
    )

    for info_field in info_fields:
        data[info_field] = list()

    for record in vcf_reader:
        print(record)
        for alt in record.ALT:
            print(alt)
            data['chromosome'].append(record.CHROM)
            data['position'].append(record.POS)
            data['reference'].append(record.REF)
            data['alternate'].append(alt)
            for info_field in info_fields:
                field_data = record.INFO[info_field]
                if info_field == 'ANTIBIOTICS':
                    result = []
                    for antibiotic_str in field_data:
                        antibiotic, obs_confidence = antibiotic_str.split('|')
                        if int(confidence) == int(obs_confidence):
                            result.append(antibiotic)
                    field_data = ','.join(result)

                data[info_field].append(field_data)

    df = pd.DataFrame.from_dict(data)
    df = df[df.ANTIBIOTICS != '']
    return df


def convert_who_confidence(who_confidence: str) -> int:
    """
    Convert WHO confidence rating.

    The WHO catalogue has a confidence rating for AMR, convert it to something
    more useful than a string with lots of characters and spaces.
    """
    confidence = {
        "1) Assoc w R": 1,
        "2) Assoc w R - Interim": 2,
        "3) Uncertain significance": 3,
        "4) Not assoc w R - Interim": 4,
        "5) Not assoc w R": 5,
        "NA": 7,
        ".": 7,
        "Synonymous": 6
    }
    return confidence[who_confidence]
