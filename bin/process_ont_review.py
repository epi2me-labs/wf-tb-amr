#!/usr/bin/env python
"""Process ONT reviewed variants."""
import Bio
from create_variant_db import (
    get_genbank_feature,
    get_reference_codon_seq,
    process_genbank,
    reverse_translate,
    translate_triplet)
import pandas as pd

review = pd.read_csv("ont_review.txt", header=0, sep="\t")

genbank_file = "../data/primer_schemes/V3/NC_000962.3.gb"

genbank = process_genbank(genbank_file)

locus_tags = {
    'Rv0678': 'Rv0678',
    'rrl': 'Rvnr02'}

print("""##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO   FORMAT  SAMPLE""")

for index, row in review.iterrows():
    locus_tag = locus_tags[row[3]]
    if row[3] == "Rv0678":
        bdq = 8
        lzd = 7
    elif row[3] == "rrl":
        bdq = 7
        lzd = 8

    codon_number = row[6]
    ref_aa = row[4]
    alt_aa = row[5]
    reference_codon = get_reference_codon_seq(genbank, locus_tag, codon_number)

    ref_codon = list()
    for i in reference_codon:
        ref_codon.append(i[1])
    ref_codon = "".join(ref_codon)

    if reverse_translate(alt_aa):
        for mutant_codon in reverse_translate(alt_aa):

            ref_three = Bio.Data.IUPACData.protein_letters_1to3[ref_aa]
            if alt_aa != '*':
                alt_three = Bio.Data.IUPACData.protein_letters_1to3[alt_aa]
                HGVS_PROTEIN = f"p.{ref_three}{codon_number}{alt_three}"
            else:
                HGVS_PROTEIN = f"p.{ref_three}{codon_number}*"
            info = dict(
                CODON_NUMBER=codon_number,
                ANTIBIOTICS=",".join([
                    'RIF|7',
                    'INH|7',
                    'EMB|7',
                    'PZA|7',
                    'LEV|7',
                    'MXF|7',
                    f"BDQ|{bdq}",
                    f"LZD|{lzd}",
                    'CFZ|7',
                    'DLM|7',
                    'AMI|7',
                    'STM|7',
                    'ETH|7',
                    'KAN|7',
                    'CAP|7']),
                GENE=row[3],
                GENE_LOCUS=locus_tag,
                FEATURE_TYPE='CDS',
                ORIGIN='WHO_CANONICAL',
                HGVS_PROTEIN=HGVS_PROTEIN,
                HGVS_NUCLEOTIDE='.',
                STRAND=[
                    feature for feature in get_genbank_feature(
                        genbank, locus_tag)
                    if feature.type == 'CDS' or feature.type == 'rRNA'
                ][0].strand
            )
            line = [
                'NC_000962.3',
                reference_codon[0][0],
                '.',
                ref_codon,
                mutant_codon,
                '.',
                '.',
                ';'.join([f'{k}={v}' for k, v in info.items()]),
                '.',
                '.'
            ]
            print("\t".join(line))

    else:
        ref_allele = row[1].upper()
        mut_allele = row[2].upper()

        count = 0
        for i in reference_codon:
            if i[1] == ref_allele:
                count += 1
        if count == 1:
            mutant_codon = []
            for i in reference_codon:
                if i[1] == row[1].upper():
                    mutant_codon.append(row[2].upper().rstrip())
                else:
                    mutant_codon.append(i[1])

            info = dict(
                CODON_NUMBER=codon_number,
                ANTIBIOTICS=",".join([
                    'RIF|7',
                    'INH|7',
                    'EMB|7',
                    'PZA|7',
                    'LEV|7',
                    'MXF|7',
                    f"BDQ|{bdq}",
                    f"LZD|{lzd}",
                    'CFZ|7',
                    'DLM|7',
                    'AMI|7',
                    'STM|7',
                    'ETH|7',
                    'KAN|7',
                    'CAP|7']),
                GENE=row[3],
                GENE_LOCUS=locus_tag,
                FEATURE_TYPE='CDS',
                ORIGIN='WHO_CANONICAL',
                HGVS_PROTEIN='.',
                HGVS_NUCLEOTIDE='.',
                STRAND=[
                    feature for feature in get_genbank_feature(
                        genbank, locus_tag)
                    if feature.type == 'CDS' or feature.type == 'rRNA'
                ][0].strand
            )
            line = [
                'NC_000962.3',
                reference_codon[0][0],
                '.',
                ref_codon,
                ''.join(mutant_codon),
                '.',
                '.',
                ';'.join([f'{k}={v}' for k, v in info.items()]),
                '.',
                '.'
            ]
            print("\t".join(line))
        else:
            features = get_genbank_feature(genbank, locus_tag)

            ref_aa = translate_triplet(ref_codon)
            mutant_codon = []
            for i in reference_codon:
                if int(features[0].location.start + row[7]) == int(i[0]):
                    if i[1] == row[1].upper():
                        mutant_codon.append(row[2].upper().rstrip())
                    else:
                        mutant_codon.append(i[1])
                else:
                    mutant_codon.append(i[1])

            mut_aa = translate_triplet(''.join(mutant_codon))
            ref3 = Bio.Data.IUPACData.protein_letters_1to3[ref_aa]
            if mut_aa != '*':
                mut3 = Bio.Data.IUPACData.protein_letters_1to3[mut_aa]
                HGVS_PROTEIN = f"p.{ref3}{codon_number}{mut3}"
            else:
                HGVS_PROTEIN = f"p.{ref3}{codon_number}*"
            info = dict(
                CODON_NUMBER=codon_number,
                ANTIBIOTICS=",".join([
                    'RIF|7',
                    'INH|7',
                    'EMB|7',
                    'PZA|7',
                    'LEV|7',
                    'MXF|7',
                    f"BDQ|{bdq}",
                    f"LZD|{lzd}",
                    'CFZ|7',
                    'DLM|7',
                    'AMI|7',
                    'STM|7',
                    'ETH|7',
                    'KAN|7',
                    'CAP|7']),
                GENE=row[3],
                GENE_LOCUS=locus_tag,
                FEATURE_TYPE='CDS',
                ORIGIN='WHO_CANONICAL',
                HGVS_PROTEIN=HGVS_PROTEIN,
                HGVS_NUCLEOTIDE='.',
                STRAND=[
                    feature for feature in get_genbank_feature(
                        genbank, locus_tag)
                    if feature.type == 'CDS' or feature.type == 'rRNA'
                ][0].strand
            )
            line = [
                'NC_000962.3',
                int(features[0].location.start + row[7]),
                '.',
                row[1].upper(),
                row[2].upper(),
                '.',
                '.',
                ';'.join([f'{k}={v}' for k, v in info.items()]),
                '.',
                '.'
            ]
            print("\t".join(line))
