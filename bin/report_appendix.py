#!/usr/bin/env python
"""Create a research report."""

import argparse
import json

from aplanat import report
from aplanat.util import ond_colors, ont_colors
from common_methods import variants_table_from_vcf


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


def make_disclaimer_section(report_doc, data):
    """Make a disclaimer section."""
    section = report_doc.add_section()
    section._add_item(f"""{make_header(data["title"])}<small><ol>""")

    for disclaimer in data['disclaimers']:
        section._add_item(f"<li>{disclaimer}</li>""")

    section._add_item("""</ol></small>""")


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
        help="WHO Variant group", default=[2, 3], nargs="+")
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
        "Mycobacterium tuberculosis appendix: further findings",
        "wf-tb-amr",
        revision=args.revision,
        commit=args.commit,
        about=False,
        style=args.style)

    with open(args.canned_text) as json_data:
        canned_text = json.load(json_data)

    section = report_doc.add_section()

    section._add_item(f"""
        <div class="mt-4">
            <p class="text-brand-primary">
                {canned_text['sections']['disclaimer']['blurb']}
            </p>
        </div>""")

    info_fields = [
        'HGVS_NUCLEOTIDE',
        'HGVS_PROTEIN',
        'GENE',
        'ANTIBIOTICS',
        'AF',
        'FS_SB']

    section._add_item(f"""
        <div class="card bg-light mt-3 mb-5">
            <h4 class="card-header">
                <i class="fa fa-clipboard-list"></i>
                Sample Details
            </h4>
            <div class="card-body">{args.sample_id}</div></div>""")

    section._add_item(f"""
        <div class="card bg-light mt-3 mb-5">
            <h4 class="card-header">
                <i class="fa fa-clipboard-list"></i>
                Tentative Additional Antibiotic Resistance
            </h4>
            <div class="card-body">
            {canned_text['appendix']['intro']}""")

    for res in args.group:
        table = variants_table_from_vcf(args.vcf_file, info_fields, res)

        section._add_item(f"""
            <div class="card bg-light mt-3">
                <h5 class="card-header">
                    <i class="fa fa-clipboard-list"></i>
                    WHO Group {res} Confidence Grade
                </h5>
                <div class="card-body">""")

        section._add_item(f"""{canned_text['who_groups'][str(res)]['name']}""")
        section._add_item("""<ol>""")
        criteria = canned_text['who_groups'][str(res)]['criteria']
        for item in criteria:
            section._add_item(f"""<li>{criteria[item]}</li>""")
        section._add_item("""</ol>""")

        section.table(
            table, index=False, key=f"group{res}",
            th_color=colors.BRAND_BLUE, paging=False, searchable=False)

        section._add_item("""</div></div>""")

    section._add_item("""</div></div>""")

    # make_disclaimer_section(
    #     report_doc, canned_text['sections']['disclaimer'])

    section = report_doc.add_section()
    section._add_item(f"""<div style="clear: both;">
        <span class="sticker" style="float: left;">RUO</span>
        <span class="sticker" style="float: left;">PROTOTYPE</span>
        <img style="float: right;" src="{colors.BRAND_LOGO}">
</div>""")

    report_doc.write(f'{args.sample_id}_appendix.html')


if __name__ == "__main__":
    main()
