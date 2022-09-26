#!/usr/bin/env python
"""Create report file."""

import argparse
import json
import math
import os
import re

from aplanat import bars, lines, report
from aplanat.components import simple as scomponents
from aplanat.util import Colors, ond_colors, ont_colors
from bokeh.layouts import gridplot, layout
from bokeh.models import (
    CategoricalColorMapper, ColumnDataSource, Legend)
from bokeh.plotting import figure
from bokeh.transform import transform
from common_methods import (
    process_sample_amr,
    process_sample_coverage,
    variants_table_from_vcf)
import numpy as np
import pandas as pd
import pyranges as pr
import pysam


def controls(
    samples: list,
    ntc_threshold: str,
    positive_threshold: str,
    report_doc
):
    """Deal with controls."""
    section = report_doc.add_section()
    section._add_item("""
    <div class="card bg-light">
        <h4 class="card-header">
            <i class="fa fa-balance-scale"></i> Run Controls
        </h4>
        <div class="card-body">
            <div class="card-group">
    """)

    no_template_stati = check_no_template(
            samples,
            ntc_threshold,
            section)

    positive_stati = check_positive(samples, positive_threshold, section)

    if no_template_stati and positive_stati:
        status = (
            "danger" if "fail" in no_template_stati
            or "fail" in positive_stati else "success")
        icon = (
            "times-circle" if "fail" in no_template_stati
            or "fail" in positive_stati else "check-circle")
        phrase = (
            "Run Failed" if "fail" in no_template_stati
            or "fail" in positive_stati else "Run Passed")
    else:
        status = "warning"
        icon = "times-circle"
        phrase = "Run status cannot be determined without proper controls."

    section._add_item(f"""
        </div>
            </div>
                <div class="card-footer">
                    <div class="alert mb-0 alert-{status}" role="alert">
                        <i class="fa fa-{icon}"></i>
                        <strong>{phrase}</strong>
                    </div>
                </div>
            </div>
    """)


def check_positive(samples: dict, threshold: str, section):
    """Output no template controls section.

    Args:
    samples (list):     list of sample dictionaries:
    {sample, type, readcount file}
    threshold (str):    threshold for failure e.g. 20,3 - where 20 is
    the number of reads and 3 is the number of
    amplicons > 20 reads would have to appear
    in to constitute a fail
    report_doc ([type]): [description]
    """
    section._add_item("""<div class="card">
    <h5 class="card-header">Positive Controls</h5>
    <div class="card-body">""")

    positive_control_count = 0
    for sample in samples:
        if samples[sample]['type'] == 'positive_control':
            positive_control_count += 1

    # if there are no positive controls then report that
    if positive_control_count == 0:
        section.alert(
            title="""<i class="fa fa-times-circle"></i>
                    No Positive Control Included""",
            text="""There was no sample labeled as
                    <strong>positive_control</strong>
                    on the samplesheet.
                    <strong>Run success cannot be determined.</strong>""",
            level="warning"
        )
        section._add_item("""</div></div>""")
        return

    # if there are positive controls then list them
    qc_stati = []

    for sample in samples:
        if samples[sample]['type'] == 'positive_control':
            qc_status = samples[sample]['qc_status']
            qc_stati.append(qc_status)
            fail_amplicons = " ".join(
                f"""<span class="badge badge-white">{k}
                 <span class="badge badge-brand-red">{int(v['median'])}</span>
                 </span>"""
                for k, v in samples[sample]['fail_amplicons'].items()
            )

            section.alert(
                title=f"""<i class="fa fa-vial"></i>
                     <a href=\"#{sample}\">{sample}</a>
                     """,
                text=(
                    f"""<p>Based on set thresholds this sample Failed</p>
                    <p>Failed Targets</p>
                    <h5>{fail_amplicons}</h5>""" if qc_status == 'fail' else
                    f"""Based on set thresholds this sample Passed
                    <h5>{fail_amplicons}</h5>"""
                    ),
                level="danger" if qc_status == 'fail' else "success"
            )

    reads, amplicons = threshold.split(',')
    section._add_item(f"""</div>
        <div class="card-footer">
            <span class="small text-muted">
                Any sample labelled as positive_control on the sample sheet
                are subject to thresholds
                <strong>
                    (<{reads} reads in >{amplicons} targets = fail)
                </strong>
                which determine test validity.
                Multiple positive_control samples can be used
            .</span>
    </div></div>""")

    return qc_stati


def check_no_template(samples: list, threshold: str, section):
    """Output no template controls section.

    Args:
    samples(list):  list of sample dictionaries:
    {sample, type, readcount file}
    threshold(str): threshold for failure e.g. 20, 3 - where 20
    is the number of reads and 3 is the number of
    amplicons > 20 reads would jhave to appear in
    to constitute a fail
    report_doc([type]): [description]
    """
    section._add_item("""<div class="card">
    <h5 class="card-header">No Template Controls</h5>
    <div class="card-body">""")

    no_template_count = 0
    for sample in samples:
        if samples[sample]['type'] == 'no_template_control':
            no_template_count += 1

    # if there are no positive controls then report that
    if no_template_count == 0:
        section.alert(
            title="""<i class="fa fa-times-circle"></i>
                No Template Control Not Included""",
            text="""There was no sample labeled as
                <strong>no_template_control</strong>
                on the sample sheet.
                <strong>Contamination can not be ruled out.</strong>""",
            level="warning"
        )
        section._add_item("""</div></div>""")
        return

    # if there are positive controls then list them
    qc_stati = []
    for sample in samples:
        if samples[sample]['type'] == 'no_template_control':
            qc_status = samples[sample]['qc_status']

            fail_amplicons = " ".join(
                f"""<span class="badge badge-white">{k}
                <span class="badge badge-brand-red">{int(v['median'])}</span>
                </span>
                """ for k, v in samples[sample]['fail_amplicons'].items())

            qc_stati.append(qc_status)

            section.alert(
                title=f"""<i class="fa fa-vial"></i>
                    <a href="#{sample}">{sample}</a>""",
                text=f"""Based on set thresholds this sample Failed
                <h5>{fail_amplicons}</h5>
                    """ if qc_status == 'fail' else f"""
                    Based on set thresholds this sample Passed
                    <h5>{fail_amplicons}</h5>""",
                level="danger" if qc_status == 'fail' else "success")

    reads, amplicons = threshold.split(',')
    section._add_item(f"""</div>
        <div class="card-footer">
            <span class="small text-muted">
                Any sample labelled as no_template_control on the
                sample sheet are subject to thresholds
                    <strong>
                        (>={reads} reads in >{amplicons} targets = fail)
                    </strong> which determine test validity.
                Multiple no_template_control samples can be used.
            </span>
        </div>
    </div>""")

    return qc_stati


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    dfs = pd.concat(dfs)
    dfs = dfs.astype({"sample_name": str})
    return dfs


def section_reads_per_barcode(args, report_doc):
    """Return a plot of reads per barcode."""
    seq_summary = read_files(args.per_barcode_stats)

    section = report_doc.add_section()

    section._add_item("""
    <div class="card bg-light mt-3">
        <h4 class="card-header">
            <i class="fa fa-poll"></i> Number of Reads Per Sample
        </h4>
        <div class="card-body">
    """)

    barcode_counts = (
        pd.DataFrame(seq_summary['sample_name'].value_counts())
        .sort_index()
        .reset_index()
        .rename(
            columns={'index': 'sample', 'sample_name': 'count'})
    )

    bc_counts = bars.simple_bar(
        barcode_counts['sample'].astype(str),
        barcode_counts['count'],
        colors=[colors.BRAND_BLUE]*len(barcode_counts),
        title=(
            'Number of reads per sample'),
        plot_width=None
    )
    bc_counts.xaxis.major_label_orientation = 3.14/2
    bc_counts.background_fill_alpha = 0
    bc_counts.border_fill_color = None
    section.plot(
        layout(
            [[bc_counts]],
            sizing_mode="stretch_width"))
    section._add_item("""</div></div>""")


def csv_output(amr, coverage, canned_text, csv):
    """Make a csv results file."""
    output = [",".join([
        "sample",
        "type",
        "status",
        "call",
        "no_resistance_detected",
        "undetermined",
        "resistant",
        "hgvs_nucleotide",
        "hgvs_protein",
        "passed_targets",
        "failed_targets"])]

    control_status = []
    for sample, items in coverage.items():
        if items['type'] in ['no_template_control', 'positive_control']:
            if items['qc_status'] == "fail":
                control_status.append("fail")

    for sample in amr:

        if coverage[sample]['qc_status'] == "pass":
            if 'fail' in control_status:
                coverage[sample]['qc_status'] = "fail"

        info = amr[sample]
        if info['type'] != 'test_sample':
            continue

        sample_type = info['type']

        line = [
            sample,
            sample_type,
            coverage[sample]['qc_status']]

        if coverage[sample]['qc_status'] != 'fail':

            abs = ";".join(list(info['resistance']['resistant'].keys()))

            nucleotides = []
            proteins = []

            for ab in list(info['resistance']['resistant'].keys()):
                ab_nuc = []
                ab_pro = []
                for variant in info['resistance']['resistant'][ab]['variants']:
                    var_info = variant['info']
                    ab_nuc.append(
                        f"{var_info['GENE']}.{var_info['HGVS_NUCLEOTIDE']}")
                    ab_pro.append(
                        f"{var_info['GENE']}.{var_info['HGVS_PROTEIN']}")

                nucleotides.append(":".join(str(x) for x in ab_nuc))
                proteins.append(":".join(str(x) for x in ab_pro))

            line = line + [
                info['resistance']['resistance_level'],
                ";".join([
                    ab for ab in info[
                        'antibiotics'] if info['antibiotics'][ab] == 0]),
                ";".join([
                    ab for ab in info[
                        'antibiotics'] if info['antibiotics'][ab] == -1]),
                abs,
                ";".join(str(x) for x in nucleotides),
                ";".join(str(x) for x in proteins)]

        else:
            line = line + ["\t" for i in range(0, 5)]

        passed_amplicons = []

        for amplicon in coverage[sample]["passed_amplicons"]:
            median = coverage[sample]['passed_amplicons'][amplicon]['median']
            passed_amplicons.append(
                f"""{amplicon}:{median}""")

        line.append(";".join(passed_amplicons))

        fail_amplicons = []

        for amplicon, fail_info in coverage[sample]["fail_amplicons"].items():
            median = fail_info['median']
            fail_amplicons.append(f"""{amplicon}:{median}""")

        line.append(";".join(fail_amplicons))

        output.append(",".join(line))

    with open(csv, "w") as f:
        f.write("\n".join(output))
    f.close()


def section_executive_summary(args, report_doc, coverage, amr, canned_text):
    """Return the results summary."""
    data = pd.DataFrame([
                {'Sample': sample, **amr[sample]['antibiotics']}
                for sample in amr if amr[sample]['type'] == 'test_sample'])

    data.set_index("Sample", inplace=True)
    data = data.replace({np.nan: 0})
    data = data.reindex(sorted(data.columns), axis=1)
    section = report_doc.add_section()
    reads, amplicons = args.sample_threshold.split(',')
    section.markdown("""
    <div class="card bg-light mt-3">
        <h4 class="card-header">
            <i class="fa fa-clipboard-list"></i> Results Summary
        </h4>
        <div class="card-body">
            <div class="card-group">
                <div class="card">
                    <h5 class="card-header">
                        <i class="fa fa-bacteria"></i> Presence/Absence of TB
                    </h5>
                    <div class="card-body">
                    <span class="badge badge-brand-primary">
                        TB present
                    </span>
                    <span class="badge badge-brand-light-grey">
                        TB present but targets missing
                    </span>
                    <span class="badge badge-brand-red">
                        TB absent or present below limit of detection
                    </span>
                    <div class="row">
                        <div class="col">""")

    # do some calculations to get total samples so we can divide into columns
    test_samples = [
        sample for sample in amr
        if amr[sample]['type'] == 'test_sample']

    total_test_samples = len(test_samples)

    half = math.ceil(total_test_samples/2)

    # keep count to know when to move on to next column
    count = 0

    for sample, info in coverage.items():

        if info['type'] != 'test_sample':
            continue

        qc_status = info['qc_status']

        # prevent failed samples from being plotted
        if qc_status == 'fail':
            data = data[data.index != sample]

        total_amplicons = info['total_amplicons']
        fail_amplicons = len(info['fail_amplicons'])
        covered_amplicons = total_amplicons-fail_amplicons
        percent_covered = (covered_amplicons / total_amplicons) * 100

        if qc_status == "pass":
            background = "success"
            progress = "brand-primary"
            symbol = '&#x2714;'
            # if there are any failed amplicons then change color to grey
            if fail_amplicons > 0:
                background = "brand-light-grey"
                progress = "brand-grey"
                symbol = '&#x203C;'

        elif qc_status == "fail":
            background = "danger"
            progress = "brand-red"
            symbol = '&#x2718;'

        section._add_item(f"""
            <div class="alert alert-{background} mt-3" role="alert">

                            {symbol} <strong><a href="#{sample}">{sample}</a>
                            <a href=\"{sample}_report.html\">
                               <i class="fa fa-clipboard"></i>
                           </a>
                            <a href=\"{sample}_appendix.html\">
                               <i class="fa fa-microscope"></i>
                           </a>
                            </strong>

                        <span class="float-right badge badge-white">
                            {amr[sample]['resistance']['resistance_level']}
                        </span>

                <div class="progress">
                    <div
                        class="progress-bar
                            progress-bar-striped
                            bg-{progress}"
                        role="progressbar"
                        style="width: {percent_covered}%;"
                        aria-valuenow="{total_amplicons-fail_amplicons}"
                        aria-valuemin="0"
                        aria-valuemax="{total_amplicons}">
                        {covered_amplicons}/{total_amplicons} Covered
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-12">
                        <small class="float-right">
                            {info['barcode']}
                        </small>
                    </div>
                </div>
        </div>""")
        count += 1

        # now we've done half the samples change to a new column
        if count == half:
            section._add_item("""</div><div class="col">""")

    section._add_item(f"""</div></div></div>
    <div class="card-footer">
        <span class="small text-muted">
            Based on set thresholds
                <strong>
                    (>={reads} reads in >{amplicons} targets)
                </strong>
            this section describes whether TB is likely present, absent or
            present below the limit of detection.
        </span>
    </div>
    </div>""")

    tools = "save"
    data.index.name = 'samples'
    data.columns.name = 'amr'
    data = data.stack().rename("value").reset_index()

    calls = {
        0: 's',  # susceptible
        -1: 'u',  # undertermined
        1: 'r1',  # resistance grp1
        2: 'r2',  # resistance grp2
        8: 'ont'}   # resistance ONT

    data["call"] = data['value'].map(calls)

    mapper = CategoricalColorMapper(
        palette=[
            colors.BRAND_LIGHT_GREY,
            '#ffffff',
            colors.BRAND_BLUE,
            aplanat_colors.sandstorm,
            aplanat_colors.fandango],
        factors=[
            's',
            'u',
            'r1',
            'r2',
            'ont'])

    p = figure(
        tools=tools,
        tooltips=None,
        plot_width=500,
        plot_height=60+(len(list(data.samples.drop_duplicates())*40)),
        y_range=list(data.samples.drop_duplicates()),
        x_range=list(data.amr.drop_duplicates()),
        toolbar_location="right",
        x_axis_location="above",
        y_axis_location="left",
        title="ONT TB Report Card"
    )
    p.rect(
        x="amr",
        y="samples",
        width=0.9,
        height=0.9,
        source=ColumnDataSource(data),
        line_color="#333333",
        fill_color=transform('call', mapper)
    )

    p.xaxis.major_label_text_font_size = "11pt"
    p.yaxis.major_label_text_font_size = "11pt"
    p.xaxis.major_label_orientation = math.pi/4
    p.title.text_font_size = "11pt"
    p.xaxis.axis_line_color = None
    p.yaxis.axis_line_color = None
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.outline_line_color = None
    p.toolbar.logo = None

    section._add_item("""
        <div class="card">
         <h5 class="card-header">
            <i class="fa fa-capsules"></i> Antimicrobial Resistance
          </h5>
          <div class="card-body">""")
    section._add_item('<div align="center"><p>')
    section.plot(p)
    section._add_item('</p></div>')
    section._add_item(f"""
         {canned_text['misc_language']['card_blurb']}
         <br><br>
         Interpretation:
            <ul>
            <li>
                <h5 style="display: inline;">
                    <span class="badge badge-brand-primary">
                        Resistant (Group 1)
                    </span>
                </h5>
            </li>
            <li>
                <h5 style="display: inline;">
                    <span class="badge"
                            style="background-color: #F5CC49; color: #ffffff;">
                        Resistant (Group 2)
                    </span>
                </h5>
            </li>
            <li>
                <h5 style="display: inline;">
                    <span class="badge"
                            style="background-color: #A53F96;color: #ffffff;">
                        Resistant (ONT)
                    </span>
                </h5>
            </li>
            <li>
                <h5 style="display: inline;">
                    <span class="badge badge-brand-light-grey">
                        No Resistance Detected
                    </span>
                </h5>
            </li>
            <li>
                <h5 style="display: inline;">
                    <span class="badge badge-white">
                        Undetermined
                    </span>
                </h5>
            </li>
        </ul>
    </div>
    </div>
    </div>
    </div>
    </div>
     """)


def target_coverage_plot_panel(
        samplekey, read_data, bed_data, ref_tiles, tile_size, report_doc,
        tab_name, ncols=4):
    """Return coverage plots."""
    plots = list()
    target_names = bed_data["tname"].to_list()
    targets = pr.PyRanges(bed_data)
    reads = pr.PyRanges(read_data)
    for i, target in enumerate(target_names):
        info = targets[targets.tname == target].df
        info.reset_index(inplace=True)

        hits = ref_tiles.count_overlaps(reads, strandedness='same')
        t = hits.intersect(pr.PyRanges(info)).df

        fwd, rev = (t[t['Strand'] == x] for x in '+-')

        names = [None, None]

        extra_width = 0 if i != ncols - 1 else 80

        p = lines.line(
            [fwd.Start.to_list(), rev.Start.to_list()],
            [fwd.NumberOverlaps, rev.NumberOverlaps+2],
            colors=[colors.BRAND_BLUE, colors.BRAND_GREY], names=names,
            height=170, width=230 + extra_width,
            ylim=(0, max(fwd.NumberOverlaps.max(),
                         rev.NumberOverlaps.max())+10),
            title='Gene: {}\n{}:{}-{}'.format(
                info['tname'][0], info['Chromosome'][0],
                info['Start'][0], info['End'][0]))
        p.xaxis.formatter.use_scientific = False
        p.xaxis.major_label_orientation = 3.14/6

        if i == ncols - 1:
            legend = Legend(
                items=[("+", p.renderers[0:1]), ("-", p.renderers[1:])])
            p.add_layout(legend, 'right')
            p.toolbar.logo = None
        plots.append(p)
    plot = gridplot(plots, ncols=ncols)

    return plot


def variants_evidence(
        args, report_doc, samples, canned_text, show_mapping_statistics=False):
    """Return variants evidence."""
    variant_canned_text = canned_text['sections']['variants']
    coverage_canned_text = canned_text['sections']['coverage']
    section = report_doc.add_section()
    section._add_item("""
        <div class="card bg-light mt-3 mb-3">
            <h5 class="card-header">
                <i class="fa fa-vials"></i> Sample Details
            </h5>
    <div class="card-body">
        The tabular information below shows a summary of
         the genetic variants observed within the sample.
         Data included in the table includes genomic coordinates, HGVS protein
         and nucleotide, the drug resistance the variant confers,
         the allele frequency of the variant (AF) and the phred scaled
         fishers exact p value of strand bias (FS_SB).
    """)
    vcf_files = args.genotype_json

    info_fields = [
        'HGVS_NUCLEOTIDE',
        'HGVS_PROTEIN',
        'GENE',
        'ANTIBIOTICS',
        'AF',
        'FS_SB']

    for vcf_file in vcf_files:

        samplekey = re.sub(
            '.final.vcf', '',
            os.path.basename(vcf_file))
        qc_status = samples[samplekey]['qc_status']
        if qc_status == 'pass':
            background = 'brand-primary'
            if len(samples[samplekey]['fail_amplicons']) > 0:
                background = "brand-light-grey"
        elif qc_status == 'fail':
            background = 'brand-red'

        section._add_item(
            f"""<div class="card mt-3 border-{background}" id={samplekey}>
            <h5 class="card-header d-flex justify-content-between
                align-items-center alert-{background}">
                <span>
                    <i class="fa fa-vial"></i>
                    {samplekey} ({samples[samplekey]['barcode']})
                </span>
                <span class="badge badge-white">
                    {samples[samplekey]['type']}
                </span>
            </h5>
            <div class="card-body">
        """)

        placeholder = report_doc.add_section(key=f"variantSummary{samplekey}")

        if qc_status == 'pass':
            section._add_item(
                f"""
                <h5>{variant_canned_text['title']}</h5>
                {variant_canned_text['blurb']}""")
            table = variants_table_from_vcf(vcf_file, info_fields, [1, 2, 8])
            placeholder.table(
                table, index=False, key=f"variantSummary{samplekey}",
                th_color=colors.BRAND_BLUE, paging=False, searchable=False)

        # recover the appropriate bed file
        bed_data = pd.read_csv(args.bed, sep="\t", header=None)
        bed_data.columns = ["Chromosome", "Start", "End", "tname"]
        bed_data["read_count"] = 0
        bed_data['region'] = [
            '{}:{}-{}'.format(x['Chromosome'], x.Start, x.End)
            for _, x in bed_data.iterrows()]

        read_data = pd.read_csv(
            f"pickedreads/{samplekey}.bamstats", sep="\t")
        read_data.rename(
            columns=dict(
                ref="Chromosome", rstart="Start",
                rend="End", direction="Strand"),
            inplace=True)

        tile_size = 20
        ref = pysam.FastaFile(args.reference)
        ref_tiles = pr.concat([
            pr.PyRanges(
                chromosomes=ref.references,
                starts=[0] * ref.nreferences,
                ends=ref.lengths,
                strands=[strand] * ref.nreferences).tile(
                    tile_size=tile_size
            ) for strand in '+-'])

        section = report_doc.add_section(
            key=f"{samplekey}-target-cover-plots")

        section._add_item(
            f"""
            <h5>{coverage_canned_text['title']}</h5>
            {coverage_canned_text['blurb']}
            {tile_size} bases.
            """)

        section.plot(target_coverage_plot_panel(
            samplekey, read_data, bed_data, ref_tiles,
            tile_size, report_doc, "Mapping data"))

        section._add_item("</div></div>")

    section._add_item("</div></div>")


def assay_details(args, report_doc, canned_text):
    """Provide some assay details."""
    section = report_doc.add_section()
    section._add_item("""
    <div class="card bg-light">
        <h4 class="card-header">
            <i class="fa-solid fa-asterisk"></i> Assay Details
        </h4>
        <div class="card-body">""")

    data = {key: list() for key in canned_text['antibiotics']['RIF'].keys()}

    data['short-name'] = list()
    for ab in canned_text['antibiotics']:
        data['short-name'].append(ab)
        for key, value in canned_text['antibiotics'][ab].items():
            if value == '':
                value = None
            data[key].append(value)

    data['genes'] = [','.join(lst) for lst in data['genes']]

    section.table(
        pd.DataFrame(data),
        th_color=colors.BRAND_BLUE, paging=False, searchable=False)

    section._add_item("</div></div>")


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
        "--per_barcode_stats", nargs='+',
        help="fastcat stats file for each sample before filtering")
    parser.add_argument(
        "--genotype_json", nargs="+",
        help="JSON file of genotypes and antibiotic resistance data")
    parser.add_argument(
        "--variants_json",
        help="JSON file of the variants under investigation")
    parser.add_argument(
        "--pickedreads", nargs="+",
        help="Mapping quality file for readgroup assigned sequences")
    parser.add_argument(
        "--reference", help="FASTA format reference genome")
    parser.add_argument(
        "--bed", help="BED file (simplified) of genomic targets")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--readcounts", default='bed_files',
        help="bedtools coverage for amplicon")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--ntc_threshold", default='20,3',
        help="comma separated string of x,y - where x is \
            the read count threshold and y is the number of amplicons \
            i.e. 20,3 - fail is more than 20 reads in more than 3 amplicons")
    parser.add_argument(
        "--sample_threshold", default='20,8',
        help="comma separated string of x,y - where x is \
            the read count threshold and y is the number of amplicons \
            i.e. 20,3 - fail is less than 20 reads in more than 8 amplicons")
    parser.add_argument(
        "--positive_threshold", default='20,2',
        help="comma separated string of x,y - where x is \
            the read count threshold and y is the number of amplicons \
            i.e. 20,2 - fail is less than 20 reads in more than 2 amplicons")
    parser.add_argument(
        "--canned_text",
        help="JSON file of text used for the report")
    parser.add_argument(
        "--versions",
        help="directory contained CSVs containing name,version.")
    parser.add_argument(
        "--style",
        help="Report style", default="ond")

    args = parser.parse_args()

    global colors
    global aplanat_colors

    if args.style == 'ond':
        colors = ond_colors
    elif args.style == 'ont':
        colors = ont_colors

    aplanat_colors = Colors

    with open(args.canned_text) as json_data:
        canned_text = json.load(json_data)

    report_doc = report.WFReport(
        "Mycobacterium tuberculosis Genotyping Run Report",
        "wf-tb-amr",
        revision=args.revision,
        commit=args.commit,
        style=args.style)

    section = report_doc.add_section()
    section._add_item(f"""
        <div class="mt-4">
            <p class="text-brand-primary">
                {canned_text['sections']['disclaimer']['blurb']}
            </p>
        </div>""")

    coverage = process_sample_coverage(
        args.metadata,
        args.readcounts,
        args.sample_threshold,
        args.ntc_threshold,
        args.positive_threshold,
        args.bed,
        canned_text)

    controls(
        coverage, args.ntc_threshold,
        args.positive_threshold, report_doc)

    amr = process_sample_amr(
        args.metadata,
        coverage,
        'variants',
        canned_text)

    csv_file = 'wf-tb-amr-report.csv'
    csv_output(amr, coverage, canned_text, csv_file)

    section_executive_summary(
        args, report_doc, coverage, amr, canned_text)

    section_reads_per_barcode(args, report_doc)

    variants_evidence(args, report_doc, coverage, canned_text)

    # lets provide some assay details
    with open(args.canned_text) as json_data:
        canned_text = json.load(json_data)
    assay_details(args, report_doc, canned_text)

    section = report_doc.add_section(
        section=scomponents.version_table(
            args.versions, th_color=colors.BRAND_BLUE))

    report_doc.add_section(
        section=scomponents.params_table(
            args.params, th_color=colors.BRAND_BLUE))

    report_doc.write('wf-tb-amr-report.html')


if __name__ == "__main__":
    main()
