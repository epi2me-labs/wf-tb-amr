#!/usr/bin/env python
"""Pytests."""

from collections import OrderedDict
import json
import tempfile

import bin.common_methods as common_methods
import bin.create_variant_db as create_variant_db
import bin.process_whatshap as process_whathap
import vcf


def test_process_resistance():
    """Test processing of variant data."""
    with open("data/primer_schemes/V2/report_config.eng.json") as json_data:
        canned_text = json.load(json_data)

    antibiotics = canned_text['antibiotics']
    vcf_file = "test_data/sample.final.variants.sorted.vcf"
    print(common_methods.process_resistance(vcf_file, antibiotics, 1))


def test_call_resistance():
    """Test calling of resistance from variant data."""
    with open("data/primer_schemes/V2/report_config.eng.json") as json_data:
        canned_text = json.load(json_data)

    antibiotics = canned_text['antibiotics']

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(),
            susceptible=OrderedDict()
        )

    resistance_result = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(
                RIF=None
            ),
            susceptible=OrderedDict()
        )

    resistance_result = OrderedDict(
            resistance_level='RR',
            resistant=OrderedDict(
                RIF=None
            ),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(
                RIF=None,
                INH=None
            ),
            susceptible=OrderedDict()
        )

    resistance_result = OrderedDict(
            resistance_level='MDR',
            resistant=OrderedDict(
                RIF=None,
                INH=None
            ),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None
            ),
            susceptible=OrderedDict()
        )

    resistance_result = OrderedDict(
            resistance_level='pre-XDR',
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None
            ),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None,
                CAP=None
            ),
            susceptible=OrderedDict()
        )

    resistance_result = OrderedDict(
            resistance_level='pre-XDR',
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None,
                CAP=None
            ),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(
                LEV=None,  # grpA fluoro
                MXF=None,  # grpA fluoro
                RIF=None,
                AMI=None,
                KAN=None,
                CAP=None,
                INH=None,
                EMB=None
            ),
            susceptible=OrderedDict()
        )

    # This one should be XDR - RIF+INH = MDR
    # MDR + 1 fluro and 1 from grpA = XDR

    resistance_result = OrderedDict(
            resistance_level='XDR',
            resistant=OrderedDict(
                LEV=None,
                MXF=None,
                RIF=None,
                AMI=None,
                KAN=None,
                CAP=None,
                INH=None,
                EMB=None
            ),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None,
                LZD=None
            ),
            susceptible=OrderedDict()
        )

    resistance_result = OrderedDict(
            resistance_level='XDR',
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None,
                LZD=None
            ),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result

    resistance = OrderedDict(
            resistance_level="NONE",
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None,
                MXF=None
            ),
            susceptible=OrderedDict()
        )

    resistance_result = OrderedDict(
            resistance_level='XDR',
            resistant=OrderedDict(
                RIF=None,
                INH=None,
                LEV=None,
                MXF=None
            ),
            susceptible=OrderedDict()
        )

    assert common_methods.call_resistance(resistance, antibiotics) == \
        resistance_result


def test_determine_status_sample_pass():
    """Test determining the status of a sample passed experiment."""
    # sample: -20,2 - fail is less than 20x in more than 2 amplicons
    sample_threshold = "-20,2"
    coverage_summary = common_methods.process_coverage(
            coverage_file="test_data/sample_pass.bedtools-coverage.bed",
            threshold=sample_threshold)

    result = {
        "status": "pass",
        "passed_targets": {
            "gene1": {
                "mean": 30.0,
                "median": 30.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 0,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 0.0,
                "passed": True
            },
            "gene2": {
                 "mean": 100.0,
                 "median": 100.0,
                 "threshold": 1.0,
                 "lt_20x_bases": 0,
                 "ge_20x_bases": 10,
                 "ge_100x_bases": 10,
                 "lt_20x_percent": 0.0,
                 "ge_20x_percent": 1.0,
                 "ge_100x_percent": 1.0,
                 "passed": True
            },
            "gene3": {
                "mean": 50.0,
                "median": 50.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 0,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 0.0,
                "passed": True
            }
        },
        "failed_targets": {
            "gene4": {
                "mean": 0.0,
                "median": 0.0,
                "threshold": 0.0,
                "lt_20x_bases": 10,
                "ge_20x_bases": 0,
                "ge_100x_bases": 0,
                "lt_20x_percent": 1.0,
                "ge_20x_percent": 0.0,
                "ge_100x_percent": 0.0,
                "passed": False
              }
           }
    }
    assert result == common_methods.determine_status(
        coverage_summary, sample_threshold)


def test_determine_status_sample_fail():
    """Test determining the status of a sample failed experiment."""
    # sample: -20,2 - fail is less than 20x in more than 2 amplicons
    sample_threshold = "-20,2"
    coverage_summary = common_methods.process_coverage(
            coverage_file="test_data/sample_fail.bedtools-coverage.bed",
            threshold=sample_threshold)

    result = {
        "status": "fail",
        "passed_targets": {
            "gene1": {
                "mean": 30.0,
                "median": 30.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 0,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 0.0,
                "passed": True
            }
        },
        "failed_targets": {
            "gene2": {
                "mean": 0.0,
                "median": 0.0,
                "threshold": 0.0,
                "lt_20x_bases": 10,
                "ge_20x_bases": 0,
                "ge_100x_bases": 0,
                "lt_20x_percent": 1.0,
                "ge_20x_percent": 0.0,
                "ge_100x_percent": 0.0,
                "passed": False
            },
            "gene3": {
                "mean": 0.0,
                "median": 0.0,
                "threshold": 0.0,
                "lt_20x_bases": 10,
                "ge_20x_bases": 0,
                "ge_100x_bases": 0,
                "lt_20x_percent": 1.0,
                "ge_20x_percent": 0.0,
                "ge_100x_percent": 0.0,
                "passed": False
            },
            "gene4": {
                "mean": 0.0,
                "median": 0.0,
                "threshold": 0.0,
                "lt_20x_bases": 10,
                "ge_20x_bases": 0,
                "ge_100x_bases": 0,
                "lt_20x_percent": 1.0,
                "ge_20x_percent": 0.0,
                "ge_100x_percent": 0.0,
                "passed": False
                }
            }
        }
    assert result == common_methods.determine_status(
        coverage_summary, sample_threshold)


def test_determine_status_ntc_fail():
    """Test determining the status of a ntc failed experiment."""
    # ntc: 20,2 - fail is more than 20x in more than 2 amplicons
    ntc_threshold = "20,2"
    coverage_summary = common_methods.process_coverage(
            coverage_file="test_data/ntc_fail.bedtools-coverage.bed",
            threshold=ntc_threshold)

    result = {
        "status": "fail",
        "failed_targets": {
            "gene1": {
                "mean": 30.0,
                "median": 30.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 0,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 0.0,
                "passed": False
            },
            "gene2": {
                "mean": 100.0,
                "median": 100.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 10,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 1.0,
                "passed": False
            },
            "gene3": {
                "mean": 50.0,
                "median": 50.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 0,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 0.0,
                "passed": False
            }
        },
        "passed_targets": {
            "gene4": {
                "mean": 0.0,
                "median": 0.0,
                "threshold": 0.0,
                "lt_20x_bases": 10,
                "ge_20x_bases": 0,
                "ge_100x_bases": 0,
                "lt_20x_percent": 1.0,
                "ge_20x_percent": 0.0,
                "ge_100x_percent": 0.0,
                "passed": True
            }
        }
    }
    assert result == common_methods.determine_status(
        coverage_summary, ntc_threshold)


def test_determine_status_ntc_pass():
    """Test determining the status of a ntc failed experiment."""
    # ntc: 20,2 - fail is more than 20x in more than 2 amplicons
    ntc_threshold = "20,2"
    coverage_summary = common_methods.process_coverage(
            coverage_file="test_data/ntc_pass.bedtools-coverage.bed",
            threshold=ntc_threshold)

    result = result = {
        "status": "pass",
        "failed_targets": {
            "gene1": {
                "mean": 30.0,
                "median": 30.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 0,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 0.0,
                "passed": False
            },
            "gene2": {
                "mean": 20.0,
                "median": 20.0,
                "threshold": 1.0,
                "lt_20x_bases": 0,
                "ge_20x_bases": 10,
                "ge_100x_bases": 0,
                "lt_20x_percent": 0.0,
                "ge_20x_percent": 1.0,
                "ge_100x_percent": 0.0,
                "passed": False
            }
        },
        "passed_targets": {
            "gene3": {
                "mean": 0.0,
                "median": 0.0,
                "threshold": 0.0,
                "lt_20x_bases": 10,
                "ge_20x_bases": 0,
                "ge_100x_bases": 0,
                "lt_20x_percent": 1.0,
                "ge_20x_percent": 0.0,
                "ge_100x_percent": 0.0,
                "passed": True
            },
            "gene4": {
                "mean": 0.0,
                "median": 0.0,
                "threshold": 0.0,
                "lt_20x_bases": 10,
                "ge_20x_bases": 0,
                "ge_100x_bases": 0,
                "lt_20x_percent": 1.0,
                "ge_20x_percent": 0.0,
                "ge_100x_percent": 0.0,
                "passed": True
            }
        }
    }
    assert result == common_methods.determine_status(
        coverage_summary, ntc_threshold)


def test_load_bed_file():
    """Test loading a bed file."""
    count = 0
    bed_file = "data/primer_schemes/V3/TB_amplicons.bed"
    with open(bed_file, "rb") as fp:
        for line in fp.readlines():
            count = count+1
    df = create_variant_db.load_bed_file(bed_file)

    assert len(df.index) == count


def test_check_position_in_bed():
    """Test check position in bed file."""
    bed_file = "data/primer_schemes/V3/TB_amplicons.bed"
    df = create_variant_db.load_bed_file(bed_file)

    pos = 7302
    check = create_variant_db.check_poistion_in_bed(pos, df)
    assert check is False

    pos = 8347
    check = create_variant_db.check_poistion_in_bed(pos, df)
    assert check is False

    pos = 8346
    check = create_variant_db.check_poistion_in_bed(pos, df)
    assert check is True

    pos = 7303
    check = create_variant_db.check_poistion_in_bed(pos, df)
    assert check is True


def test_convert_who_confidence():
    """Test convert WHO confidence ratings."""
    assert create_variant_db.convert_who_confidence("5) Not assoc w R") == 5
    assert create_variant_db.convert_who_confidence("1) Assoc w R") == 1


def test_get_reference_codon_seq_positive():
    """Test getting reference codon from a gene on the positive strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    # gyrB_L375V
    locus_tag = 'Rv0005'
    codon_number = 375
    result = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)
    assert result == [(6362, 'T'), (6363, 'T'), (6364, 'G')]


def test_get_reference_codon_seq_negative():
    """Test getting reference codon from a gene on the negative strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    # gyrB_L375V
    locus_tag = 'Rv3919c'
    codon_number = 192
    result = create_variant_db.get_reference_codon_seq(
        genbank,
        locus_tag,
        codon_number
    )
    print(result)
    assert result == [(4407629, 'G'), (4407628, 'G'), (4407627, 'C')]


def test_make_mutant_codon_positive():
    """Test making a mutant codon from a gene on the positive strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    locus_tag = 'Rv0005'
    codon_number = 375
    reference_codon = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)
    position = 6364
    alt_base = 'C'
    result = create_variant_db.make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )
    assert result == [(6362, 'T'), (6363, 'T'), (6364, 'C')]


def test_make_mutant_codon_positive_multiple():
    """Test making a mutant codon from a gene on the positive strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    locus_tag = 'Rv0005'
    codon_number = 375
    reference_codon = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)
    position = 6363
    alt_base = 'GG'
    result = create_variant_db.make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )
    print(result)
    assert result == [(6362, 'T'), (6363, 'G'), (6364, 'G')]


def test_make_mutant_codon_positive_multiple2():
    """Test making a mutant codon from a gene on the positive strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    locus_tag = 'Rv0005'
    codon_number = 375
    reference_codon = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)
    position = 6362
    alt_base = 'AAA'
    result = create_variant_db.make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )
    print(result)
    assert result == [(6362, 'A'), (6363, 'A'), (6364, 'A')]


def test_make_mutant_codon_negative():
    """Test making a mutant codon from a gene on the negative strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    locus_tag = 'Rv3919c'
    codon_number = 192
    reference_codon = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)
    position = 4407629
    alt_base = 'G'
    result = create_variant_db.make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )
    assert result == [(4407629, 'C'), (4407628, 'G'), (4407627, 'C')]


def test_make_mutant_codon_negative_multiple():
    """Test making a mutant codon from a gene on the negative strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    locus_tag = 'Rv0676c'
    codon_number = 948
    reference_codon = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)
    print(reference_codon)
    position = 775638
    alt_base = 'GC'
    result = create_variant_db.make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )
    print(result)

    assert result == [(775639, 'G'), (775638, 'C'), (775637, 'T')]


def test_make_mutant_codon_negative_multiple2():
    """Test making a mutant codon from a gene on the negative strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    locus_tag = 'Rv3919c'
    codon_number = 192
    reference_codon = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)
    position = 4407627
    alt_base = 'AAT'

    result = create_variant_db.make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )
    assert result == [(4407629, 'A'), (4407628, 'T'), (4407627, 'T')]


def test_make_mutant_codon_ubiA():
    """Test making a mutant codon from a gene on the negative strand."""
    genbank = create_variant_db.process_genbank(
        "data/primer_schemes/V2/NC_000962.3.gb")
    locus_tag = 'Rv3806c'
    codon_number = 249
    reference_codon = create_variant_db.get_reference_codon_seq(
        genbank, locus_tag, codon_number)

    position = 4269088
    alt_base = 'A'
    result = create_variant_db.make_mutant_codon(
        genbank,
        locus_tag,
        reference_codon,
        position,
        alt_base
    )
    print(result)

    assert result == [(4269089, 'G'), (4269088, 'T'), (4269087, 'G')]


def test_variants_table_from_vcf():
    """Test getting variants from VCF."""
    info_fields = [
        'HGVS_NUCLEOTIDE',
        'HGVS_PROTEIN',
        'GENE',
        'ANTIBIOTICS',
        'AF',
        'STRAND_BIAS']

    print(common_methods.variants_table_from_vcf(
        "test_data/sample.final.variants.sorted.vcf", info_fields, 1))


def test_process_whatshap():
    """Test processing whatshap results."""
    tmp = tempfile.NamedTemporaryFile(delete=False)

    phased_vcf = "test_data/phasing_test_input.vcf"
    truth_vcf = "test_data/phasing_test_output.vcf"
    template = "data/template.vcf"
    processed_vcf = tmp.name

    process_whathap.process_whathap(phased_vcf, template, processed_vcf)

    vcf_truth = vcf.Reader(filename=truth_vcf)

    vcf_test = vcf.Reader(filename=processed_vcf)

    truth = dict()
    test = dict()

    for record in vcf_test:
        test[record.POS] = record

    for record in vcf_truth:
        truth[record.POS] = record

    for position in truth:
        assert truth[position].REF == test[position].REF
        assert truth[position].ALT == test[position].ALT

    for position in test:
        assert truth[position].REF == test[position].REF
        assert truth[position].ALT == test[position].ALT
