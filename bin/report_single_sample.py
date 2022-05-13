#!/usr/bin/env python
"""Create single sample report."""

import argparse
import base64
import datetime
import io
import json

from aplanat import report
from aplanat.util import ond_colors, ont_colors
from common_methods import (
    call_resistance, comma_separator, determine_status,
    process_coverage, process_resistance)
import qrcode


def make_header(section_name: str):
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
        sample_id: str, barcode: str,
        section: report.HTMLSection, resistance: dict):
    """Make a QR code from resistance data."""
    qr = qrcode.QRCode(
        version=4,
        error_correction=qrcode.constants.ERROR_CORRECT_H,
        box_size=2,
        border=0)

    res = '|'.join([antibiotic for antibiotic in resistance['resistant']])

    level = resistance['resistance_level']
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
        report_doc: report, data: dict, sample_id: str,
        resistance: dict, barcode: str):
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
        report_doc: report, data: dict, barcode: str,
        revision: str, commit: str):
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
        report_doc: report, coverage_file: str, ntc_coverage_file: str,
        pos_coverage_file: str, ntc_threshold: str, pos_threshold: str,
        sample_threshold: str, data: dict):
    """Assess coverage of controls and sample."""
    tests = list()
    section = report_doc.add_section()
    section._add_item(f"""
        {make_header(data["title"])}
        <div class="row">""")

    controls = dict(
        sample=dict(threshold=sample_threshold, coverage_file=coverage_file),
        ntc=dict(threshold=ntc_threshold, coverage_file=ntc_coverage_file),
        pos=dict(threshold=pos_threshold, coverage_file=pos_coverage_file))

    for k, v in controls.items():

        coverage = process_coverage(v['coverage_file'], v['threshold'])

        result = determine_status(coverage, v['threshold'])

        # failed_targets = [
        #     f"""<span class="badge badge-danger">{target}
        #             <span class="badge badge-light">
        #                 {int(result['failed_targets'][target]['median'])}
        #             </span>
        #         </span>""" for target in result['failed_targets']]

        if result['status'] == 'pass':

            section._add_item(f"""
                <div class="col-4">
                    <div class="alert alert-brand-primary">
                        &#x2714; {k.upper()} status: PASS
                    </div>
                </div>""")

            tests.append(True)

        if result['status'] == 'fail':

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

    return True


def make_final_result_section(
        section: report.HTMLSection, data: dict, antibiotics: dict,
        resistance: dict, sample_id: str):
    """Make a final result section, a text describing the result."""
    antibiotics_data = list()

    for antibiotic in resistance['resistant']:
        if len(resistance['resistant'][antibiotic]['variants']) == 0:
            continue

        antibiotics_data.append(
            antibiotics[antibiotic]['full-name'])

    final_text = data["result"]["positive"].format(
        resistance=comma_separator(antibiotics_data), sample_id=sample_id)

    section._add_item(f"""
        <div class="alert alert-brand-primary">{final_text}</div>""")


def make_lineage_section():
    """Make a lineage section - we don't do this yet."""
    pass


def make_drug_section(
        report_doc: report, data: dict, antibiotics: dict, resistance: dict,
        sample_id: str, barcode: str):
    """Make the detailed drug/variant section."""
    section = report_doc.add_section()

    section._add_item(
        f"""{make_header(data['sections']['drug_susceptibility']["title"])}""")

    make_final_result_section(
        section, data['sections']['final'], antibiotics, resistance, sample_id)

    section._add_item(f"""
        <div class="row">
            <div class="col">
                {data['sections']['drug_susceptibility']['blurb']}
            </div>
        <div class="col">""")

    for k, v in data['sections']['drug_susceptibility']['result'].items():
        checked = ''
        if k == resistance["resistance_level"]:
            checked = " checked"
            v = f"<strong>{v}</strong>"

        section._add_item(f"""
            <input type=\"checkbox\" style=\"pointer-events: none;\"{checked}>
            {v}<br>""")

    section._add_item("""</div>""")

    qr_string, b64 = make_result_qr_code(
        sample_id, barcode, section, resistance)

    section._add_item(
        f"""<div class="col-2">
        <img src="data:image/png;base64, {b64}" />
        </div></div>
        <div class="row">""")

    color = dict(resistant="td-brand-primary", susceptible="td-brand-grey")
    symbol = dict(resistant="&#x2731;", susceptible="")

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

        for status in ['resistant', 'susceptible']:
            interpretation = data['misc_language'][status]
            for antibiotic in resistance[status]:
                if line != resistance[status][antibiotic]['line']:
                    continue
                gene_targets = list()
                if len(resistance[status][antibiotic]['variants']) != 0:
                    for variant in resistance[status][antibiotic]['variants']:
                        print(variant.INFO['HGVS_PROTEIN'])
                        hgvs = (
                            variant.INFO['HGVS_PROTEIN']
                            if variant.INFO['HGVS_PROTEIN'] is not None
                            else variant.INFO['HGVS_NUCLEOTIDE'])
                        gene_targets.append(
                            f"""{variant.INFO['GENE']} ({hgvs},
                            {int(variant.INFO['AF']*100)}%)""")

                section._add_item(f"""
                    <tr class=\"{color[status]}\">
                    <td>{symbol[status]} {interpretation}</td>
                    <td>{antibiotics[antibiotic]['full-name'].capitalize()}</td>
                    <td>{'<br>'.join(gene_targets)}</td>
                    </tr>""")

        section._add_item("""</table></div>""")

    section._add_item("""</div>""")


def make_disclaimer_section(report_doc: report, data: dict):
    """Make a disclaimer section."""
    section = report_doc.add_section()
    section._add_item(f"""{make_header(data["title"])}<small><ol>""")

    for disclaimer in data['disclaimers']:
        section._add_item(f"<li>{disclaimer}</li>""")

    section._add_item("""</ol></small>""")


def make_authorisation_section(report_doc: report, data: dict):
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
        "--sample_coverage_file",
        help="Barcode used for sample")
    parser.add_argument(
        "--ntc_coverage_file",
        help="Barcode used for sample")
    parser.add_argument(
        "--pos_coverage_file",
        help="Barcode used for sample")
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

    resistance = process_resistance(
        args.vcf_file, canned_text['antibiotics'], args.group)

    resistance = call_resistance(resistance, canned_text['antibiotics'])

    make_sample_section(
        report_doc, canned_text['sections']['sample'],
        args.sample_id, resistance, args.barcode)

    make_assay_section(
        report_doc, canned_text['sections']['assay'], args.barcode,
        args.revision, args.commit)

    controls = make_controls_section(
        report_doc, args.sample_coverage_file, args.ntc_coverage_file,
        args.pos_coverage_file, args.ntc_threshold, args.pos_threshold,
        args.sample_threshold, canned_text['sections']['controls'])

    # something has has failed coverage checks - make report and exit
    if controls is True:

        make_drug_section(
            report_doc, canned_text, canned_text['antibiotics'],
            resistance, args.sample_id, args.barcode)

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
