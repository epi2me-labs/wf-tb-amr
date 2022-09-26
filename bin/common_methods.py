#!/usr/bin/env python
"""Methods common to reporting."""
import copy
import json

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
    print(res)
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
        vcf_file: dict, antibiotics: dict, set_scores: list) -> dict:
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
            score = int(score)
            if score in set_scores:
                if antibiotic not in result['resistant']:
                    result['resistant'][antibiotic] = \
                        copy.deepcopy(antibiotics[antibiotic])
                    result['resistant'][antibiotic]['groups'] = list()
                    result['resistant'][antibiotic]['variants'] = list()

                if score not in result['resistant'][antibiotic]['groups']:
                    result['resistant'][antibiotic]['groups'].append(score)

                result['resistant'][antibiotic]['variants'].append(dict(
                    pos=record.POS,
                    ref=record.REF,
                    alt=str(record.ALT[0]),
                    info=dict(record.INFO)))

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
        coverage_file: str, threshold: str, bed: str) -> dict:
    """Process coverage data for a sample."""
    try:
        coverage = pd.read_csv(coverage_file, sep='\t', header=None)
    except FileNotFoundError:
        # it is expected sometimes that a NTC may have no reads
        # in this case we will create a summary with 0
        amplicons = pd.read_csv(bed, sep='\t', header=None)
        summary = dict()
        for gene in amplicons[3]:
            summary[gene] = dict(
                mean=0,
                median=0,
                threshold=0,
                lt_20x_bases=0,
                ge_20x_bases=0,
                ge_100x_bases=0,
                lt_20x_percent=0,
                ge_20x_percent=0,
                ge_100x_percent=0
            )
        summary = pd.DataFrame.from_dict(summary, orient='index')
        return summary

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


def determine_status(
        coverage_summary: pd.DataFrame,
        threshold: str,
        canned_text: dict) -> dict:
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

    # process failed antibiotics?
    genes = {
        gene: []
        for ab in canned_text['antibiotics']
        for gene in canned_text['antibiotics'][ab]['genes']}

    for ab in canned_text['antibiotics']:
        for gene in canned_text['antibiotics'][ab]['genes']:
            genes[gene].append(ab)

    for failed_target in result['failed_targets']:
        if failed_target in genes:
            result['failed_targets'][failed_target]['antibiotics'] = genes[
                failed_target]

    for pass_target in result['passed_targets']:
        if pass_target in genes:
            result['passed_targets'][pass_target]['antibiotics'] = genes[
                pass_target]

    return result


def variants_table_from_vcf(
        vcf_file: str, info_fields: list, confidence: list) -> pd.DataFrame:
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
        for alt in record.ALT:
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
                        if int(obs_confidence) in confidence:
                            result.append(f"{antibiotic} ({obs_confidence})")
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


def process_sample_coverage(
        metadata,
        readcounts,
        sample_threshold,
        ntc_threshold,
        positive_threshold,
        bed,
        canned_text) -> dict:
    """Process inputs and make some decisions on QC."""
    bed_ext = "bedtools-coverage.bed"

    with open(metadata) as metadata:
        sample_coverage = {
            d['sample_id']: {
                'type': d['type'],
                'barcode': d['barcode'],
                'readcount': f"{readcounts}/{d['sample_id']}.{bed_ext}"
            } for d in json.load(metadata)
        }

    controls = dict(
        test_sample=dict(threshold=sample_threshold),
        no_template_control=dict(threshold=ntc_threshold),
        positive_control=dict(threshold=positive_threshold))

    for sample in sample_coverage:

        sample_type = sample_coverage[sample]['type']

        coverage = process_coverage(
            sample_coverage[sample]['readcount'],
            controls[sample_type]['threshold'],
            bed)

        result = determine_status(
            coverage,
            controls[sample_type]['threshold'],
            canned_text)

        sample_coverage[sample]['total_covered'] = len(
            result['passed_targets'])

        sample_coverage[sample]['qc_status'] = result['status']
        sample_coverage[sample]['fail_amplicons'] = result['failed_targets']
        sample_coverage[sample]['passed_amplicons'] = result['passed_targets']
        sample_coverage[sample]['total_amplicons'] = (
            len(result['passed_targets'])+len(result['failed_targets']))

    return sample_coverage


def process_sample_amr(
        metadata: str,
        coverage: dict,
        variant_dir: str,
        canned_text: dict):
    """For a set of samples process AMR."""
    vcf_ext = 'final.vcf'

    with open(metadata) as metadata:
        sample_amr = {
            d['sample_id']: {
                'type': d['type'],
                'barcode': d['barcode'],
                'vcf': f"{variant_dir}/{d['sample_id']}.{vcf_ext}"
            } for d in json.load(metadata)
        }
    metadata.close()

    for sample in sample_amr:

        info = sample_amr[sample].copy()

        sample_type = info['type']

        if sample_type in ['no_template_control', 'positive_control']:
            with open(f"jsons/{sample}_data.json", 'w') as j:
                j.write(
                    json.dumps(
                        {sample: {'coverage': coverage[sample]}}, indent=4))
                continue

        processed_resistance = process_resistance(
            info['vcf'], canned_text['antibiotics'], [1, 2, 8])

        resistance = call_resistance(
            processed_resistance, canned_text['antibiotics'])

        sample_amr[sample]['resistance'] = resistance

        all_antibiotics = (
            [antibiotic for antibiotic in resistance['resistant']] +
            [antibiotic for antibiotic in resistance['susceptible']])

        sample_amr[sample]['antibiotics'] = dict()

        for antibiotic_name in all_antibiotics:
            sample_amr[sample]['antibiotics'][antibiotic_name] = 0

        for antibiotic in resistance['resistant']:
            sample_amr[sample]['antibiotics'][antibiotic] = min(
                resistance['resistant'][antibiotic]['groups'])

        # set those with failed coverage to undetermined
        failed_amplicons = coverage[sample]['fail_amplicons']
        for amplicon in failed_amplicons:
            if 'antibiotics' in failed_amplicons[amplicon]:
                for ab in failed_amplicons[amplicon]['antibiotics']:
                    # if resistance has been found in another amplicon
                    # then do not set ab to undetermined
                    if sample_amr[sample]['antibiotics'][ab] not in [1, 2, 8]:
                        sample_amr[sample]['antibiotics'][ab] = -1

        # write results to json for use in other downstream reporting steps
        with open(f"jsons/{sample}_data.json", 'w') as j:
            j.write(
                json.dumps(
                    {sample: {
                        'amr': sample_amr[sample],
                        'coverage': coverage[sample]}}, indent=4))

    return sample_amr
