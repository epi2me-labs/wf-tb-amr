#!/usr/bin/env python
"""Create single sample report."""

import argparse
import base64
import datetime
import io
import json

from aplanat import report
from aplanat.util import ond_colors, ont_colors
from common_methods import comma_separator
import qrcode


def make_header(section_name):
    """Make a section header."""
    header = f"""
    <div class="row no-gutters mt-3">
        <div class="col-1  pd-0">
            <hr style="border-top: 1px solid black;">
        </div>
        <div class="col-md-auto">
            <h5 class="ml-2 mr-2">{section_name}</h5>
        </div>
        <div class="col  pd-0">
            <hr style="border-top: 1px solid black;">
        </div>
    </div>"""

    return header


def make_result_qr_code(
        sample_id, barcode,
        section, resistance):
    """Make a QR code from resistance data."""
    qr = qrcode.QRCode(
        version=4,
        error_correction=qrcode.constants.ERROR_CORRECT_H,
        box_size=2,
        border=0)

    res = '|'.join([
        antibiotic for antibiotic, value in resistance[
            'antibiotics'].items() if value in [1, 2, 8]])

    level = resistance['resistance']['resistance_level']
    qr_string = f"""{sample_id},{barcode},{level},RES;{res}"""

    qr.add_data(qr_string)

    qr.make()

    img = qr.make_image()

    with io.BytesIO() as output:
        img.save(output, format="png")
        contents = output.getvalue()

    b64 = base64.b64encode(contents).decode("utf-8")

    return (qr_string, b64)


def make_sample_section(
        report_doc, data, sample_id, resistance, barcode):
    """Make a section with sample details."""
    section = report_doc.add_section()
    section._add_item(f"""
        {make_header(data["title"])}
        <div class="row">
            <div class="col-12">
                <table class="table border-0 table-sm mb-0"><thead><tr>""")

    for k, v in data['table'].items():
        section._add_item(f"""<th>{k}</th>""")

    section._add_item("""</tr><thead><tr>""")

    values = dict(
        sample_id=sample_id,
        report_date=datetime.datetime.now().strftime("%Y-%m-%d"))

    for k, v in data['table'].items():
        value = v.format(**values)
        section._add_item(f"<td>{value}</td>")

    section._add_item("""</tr></table></div></div>""")


def make_assay_section(
        report_doc, data, barcode, revision, commit):
    """Make a section with the assay details."""
    section = report_doc.add_section()
    section._add_item(f"""
        {make_header(data["title"])}
        <div class="row no-gutters">
            <div class="col">
            <table class="table table-sm mb-0">
                <thead>
                    <tr>""")

    for k, v in data['table'].items():
        section._add_item(f"""<th>{k}</th>""")

    section._add_item("""</tr><thead><tr>""")

    values = dict(
        platform="",
        version=f"{revision} ({commit})",
        reference="NC_000962.3",
        barcode=barcode
    )

    for k, v in data['table'].items():
        value = v.format(**values)
        section._add_item(f"<td>{value}</td>")

    section._add_item("""</tr></table></div></div>""")


def make_controls_section(
        sample, report_doc, coverage, canned_text):
    """Assess coverage of controls and sample."""
    canned_statements = canned_text['sections']['controls']

    tests = list()

    section = report_doc.add_section()

    section._add_item(f"""
        {make_header(canned_statements["title"])}
        <div class="row">""")

    controls = dict(
        ntc=[
            data['coverage'] for sample_id, data in coverage.items() if data[
                'coverage']['type'] == "no_template_control"][0],
        pos=[
            data['coverage'] for sample_id, data in coverage.items() if data[
                'coverage']['type'] == "positive_control"][0],
        sample=[
            data['coverage'] for sam_id, data in coverage.items() if data[
                'coverage']['type'] == "test_sample" and sam_id == sample][0]
    )

    for k, v in controls.items():

        if v['qc_status'] == 'pass':

            section._add_item(f"""
                <div class="col-4">
                    <div class="alert alert-brand-primary">
                        &#x2714; {k.upper()} status: PASS
                    </div>
                </div>""")

            tests.append(True)

        if v['qc_status'] == 'fail':

            section._add_item(f"""
            <div class="col-4">
                <div class="alert alert-brand-light-grey">
                    &#x2718; {k.upper()} status: FAIL
                    </div>
                </div>""")

            tests.append(False)

    section._add_item("""</div>""")

    if False in tests:
        return False

    return controls


def make_final_result_section(
        section, data, antibiotics, resistance, sample_id):
    """Make a final result section, a text describing the result."""
    antibiotics_data = list()

    resistance_resistant = resistance[
        sample_id]['amr']['resistance']['resistant']

    for antibiotic in resistance_resistant:
        antibiotics_data.append(
            antibiotics[antibiotic]['full-name'])

    final_text = data["result"]["positive"].format(
        resistance=comma_separator(antibiotics_data), sample_id=sample_id)

    section._add_item(f"""
        <div class="alert alert-brand-primary">{final_text}</div>""")

    failed_abs = [
        ab for ab, value in resistance[
            sample_id]['amr']['antibiotics'].items() if value == -1]
    if len(failed_abs) > 0:
        anbs = [antibiotics[ab]['full-name'] for ab in failed_abs]
        final_text = f"""Targets failed for {comma_separator(anbs)}"""
        section._add_item(f"""
            <div class="alert alert-brand-red">{final_text}</div>""")


def make_lineage_section():
    """Make a lineage section - we don't do this yet."""
    pass


def make_drug_section(
        report_doc, data, antibiotics, all_info,
        sample_id, barcode, coverage_result):
    """Make the detailed drug/variant section."""
    section = report_doc.add_section()

    section._add_item(
        f"""{make_header(data['sections']['drug_susceptibility']["title"])}""")

    make_final_result_section(
        section,
        data['sections']['final'],
        antibiotics,
        all_info,
        sample_id)

    section._add_item(f"""
        <div class="row">
            <div class="col">
                {data['sections']['drug_susceptibility']['blurb']}
            </div>
        <div class="col">""")

    for k, v in data['sections']['drug_susceptibility']['result'].items():
        checked = ''
        if k == all_info[sample_id]['amr']['resistance']["resistance_level"]:
            checked = " checked"
            v = f"<strong>{v}</strong>"

        section._add_item(f"""
            <input type=\"checkbox\" style=\"pointer-events: none;\"{checked}>
            {v}<br>""")

    section._add_item("""</div>""")

    qr_string, b64 = make_result_qr_code(
        sample_id, barcode, section, all_info[sample_id]['amr'])

    section._add_item(
        f"""<div class="col-2">
        <img src="data:image/png;base64, {b64}" />
        </div></div>
        <div class="row">""")

    color = dict(
        resistant="td-brand-primary",
        susceptible="td-brand-grey",
        target_fail="td-brand-red")

    symbol = dict(
        resistant="&#x2731;",
        susceptible="",
        target_fail="")

    for line in ['first-line', 'second-line']:
        section._add_item(f"""<div class="col-6">
        <table class="table table-sm mt-4 mb-1">
            <thead>
                <tr>
                    <th>Interpretation</th>
                    <th>Drug</th>
                    <th>Gene (Mutation, Allele %)</th>
                </tr>
            </thead>
            <tr>
                <td colspan=4 class="">
                    <strong>
                        {data['sections']['drug_susceptibility']["line"][line]}
                    </strong>
                </td>
        </tr>""")

        resistance = all_info[sample_id]['amr']['resistance']['resistant']
        status_codes = {
            'resistant': [1, 2, 8],
            'susceptible': [0],
            'target_fail': [-1]}
        for status, codes in status_codes.items():
            interpretation = data['misc_language'][status]
            anbs = all_info[sample_id]['amr']['antibiotics']
            for antibiotic, code in anbs.items():
                if code not in codes:
                    continue

                if line != data["antibiotics"][antibiotic]['line']:
                    continue

                gene_targets = list()
                if antibiotic in resistance:

                    for variant in resistance[antibiotic]['variants']:
                        hgvs = (
                            variant['info']['HGVS_PROTEIN']
                            if variant['info']['HGVS_PROTEIN'] is not None
                            else variant['info']['HGVS_NUCLEOTIDE'])
                        gene_targets.append(
                            f"""{variant['info']['GENE']} ({hgvs},
                            {int(variant['info']['AF']*100)}%)""")

                section._add_item(f"""
                    <tr class=\"{color[status]}\">
                    <td>{symbol[status]} {interpretation}</td>
                    <td>{antibiotics[antibiotic]['full-name'].capitalize()}</td>
                    <td>{'<br>'.join(gene_targets)}</td>
                    </tr>""")

        section._add_item("""</table></div>""")

    section._add_item("""</div>""")


def make_disclaimer_section(report_doc, data):
    """Make a disclaimer section."""
    section = report_doc.add_section()
    section._add_item(f"""{make_header(data["title"])}<small><ol>""")

    for disclaimer in data['disclaimers']:
        section._add_item(f"<li>{disclaimer}</li>""")

    section._add_item("""</ol></small>""")


def make_authorisation_section(report_doc, data):
    """Make an authorisation section."""
    section = report_doc.add_section()
    section._add_item(f"""
        {make_header(data["title"])}
        <table class="table thead-light table-sm mb-0">
            <thead>
                <tr>""")

    for k, v in data['table'].items():
        section._add_item(f"""<th>{k}</th>""")

    section._add_item("""</thead></tr><tr>""")
    for k, v in data['table'].items():
        section._add_item(f"""<td>{v}</td>""")

    section._add_item("""</tr></table>""")


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser(
        "Microbial genotyping QC report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)

    parser.add_argument(
        "--revision", default='unknown',
        help="revision number")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit number")
    parser.add_argument(
        "--canned_text",
        help="JSON file of text used for the report")
    parser.add_argument(
        "--sample_id",
        help="Sample identifier")
    parser.add_argument(
        "--barcode",
        help="Barcode used for sample")
    parser.add_argument(
        "--jsons", nargs="+",
        help="model for results")
    parser.add_argument(
        "--ntc_threshold",
        help="Threshold for ntc coverage")
    parser.add_argument(
        "--pos_threshold",
        help="Threshold for positive coverage")
    parser.add_argument(
        "--sample_threshold",
        help="Threshold for sample coverage")
    parser.add_argument(
        "--vcf",
        dest="vcf_file",
        help="VCF of the final sample variants.")
    parser.add_argument(
        "--group",
        help="WHO Variant group", default=1, type=int)
    parser.add_argument(
        "--style",
        help="Report style", default="ond")

    args = parser.parse_args()

    global colors

    if args.style == 'ond':
        colors = ond_colors
    elif args.style == 'ont':
        colors = ont_colors

    report_doc = report.WFReport(
        "Mycobacterium tuberculosis Genotyping Report",
        "wf-tb-amr",
        revision=args.revision,
        commit=args.commit,
        about=False,
        style=args.style)

    with open(args.canned_text) as json_data:
        canned_text = json.load(json_data)

    # with open(args.genotype_json) as json_data:
    #     variants_data = json.load(json_data)

    section = report_doc.add_section()

    section._add_item(f"""
        <div class="mt-4">
            <p class="text-brand-primary">
                {canned_text['sections']['disclaimer']['blurb']}
            </p>
        </div>""")

    jsons = list()
    for json_file in args.jsons:
        with open(json_file) as json_data:
            data = json.load(json_data)
            jsons.append(data)

    all_info = dict()
    for data in jsons:
        all_info = {**all_info, **data}

    make_sample_section(
        report_doc, canned_text['sections']['sample'],
        args.sample_id, all_info[args.sample_id]['amr'], args.barcode)

    make_assay_section(
        report_doc, canned_text['sections']['assay'], args.barcode,
        args.revision, args.commit)

    coverage_result = make_controls_section(
        args.sample_id, report_doc, all_info, canned_text)

    # something has has failed coverage checks - make report and exit
    if coverage_result is not False:

        make_drug_section(
            report_doc, canned_text, canned_text['antibiotics'],
            all_info, args.sample_id, args.barcode, coverage_result)

    make_authorisation_section(
        report_doc, canned_text['sections']['authorisation'])

    make_disclaimer_section(
        report_doc, canned_text['sections']['disclaimer'])

    section = report_doc.add_section()
    section._add_item(f"""<div style="clear: both;">
        <span class="sticker" style="float: left;">RUO</span>
        <span class="sticker" style="float: left;">PROTOTYPE</span>
        <img style="float: right;" src="{colors.BRAND_LOGO}">
        </div>""")

    report_doc.write(f'{args.sample_id}_report.html')


if __name__ == "__main__":
    main()
