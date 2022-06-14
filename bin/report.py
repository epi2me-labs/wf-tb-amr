#!/usr/bin/env python
"""Create report file."""

import argparse
import json
import math
import os
import re

from aplanat import bars, lines, report
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
import aplanat.graphics
from aplanat.util import ond_colors, ont_colors
from bokeh.layouts import gridplot, layout
from bokeh.models import ColumnDataSource, Legend, LinearColorMapper, Panel
from bokeh.plotting import figure
from bokeh.transform import transform
from common_methods import (
    call_resistance, determine_status,
    process_coverage, process_resistance, variants_table_from_vcf)
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


def fastcat_report_tab(file_name, tab_name):
    """Read fastcat dataframe and create a tab with qual and len plots."""
    df = pd.read_csv(file_name, sep='\t')
    depth = len(df.index)
    min_length = df["read_length"].min()
    max_length = df["read_length"].max()
    lengthplot = fastcat.read_length_plot(
        df,
        min_len=min_length,
        max_len=max_length)
    qstatplot = fastcat.read_quality_plot(df)
    exec_summary = aplanat.graphics.InfoGraphItems()
    exec_summary.append(
        'No. reads', str(depth), "bars", '')
    exec_plot = aplanat.graphics.infographic(
        exec_summary.values(), ncols=1)
    tab = Panel(
        child=layout(
            [[exec_plot], [lengthplot], [qstatplot]], aspect_ratio="auto",
            sizing_mode='stretch_width'), title=tab_name)
    return tab


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


def csv_output(samples, canned_text, csv):
    """Make a csv results file."""
    output = ["sample,type,status,call,resistant,hgvs_nucleotide,hgvs_protein"]

    for sample in samples:
        # print(samples[sample])
        sample_type = samples[sample]['type']

        vcf_file = f"variants/{sample}.final.vcf"

        resistance = process_resistance(
            vcf_file, canned_text['antibiotics'], 1)

        resistance = call_resistance(
            resistance, canned_text['antibiotics'])

        abs = ";".join(list(resistance['resistant'].keys()))
        nuc = []
        pro = []
        for ab in list(resistance['resistant'].keys()):
            ab_nuc = []
            ab_pro = []
            for variant in resistance['resistant'][ab]['variants']:
                info = variant.INFO
                ab_nuc.append(
                    f"{info['GENE']}.{info['HGVS_NUCLEOTIDE']}")
                ab_pro.append(
                    f"{info['GENE']}.{info['HGVS_PROTEIN']}")
            nuc.append(":".join(str(x) for x in ab_nuc))
            pro.append(":".join(str(x) for x in ab_pro))

        line = [
            sample,
            sample_type,
            samples[sample]['qc_status'],
            resistance['resistance_level'],
            abs,
            ";".join(str(x) for x in nuc),
            ";".join(str(x) for x in pro)
            ]

        output.append(",".join(line))

    with open(csv, "w") as f:
        f.write("\n".join(output))
    f.close()


def section_executive_summary(args, report_doc, samples, canned_text):
    """Return the results summary."""
    result = dict()

    for sample in samples:
        sample_type = samples[sample]['type']
        if sample_type in ['no_template_control', 'positive_control']:
            continue

        vcf_file = f"variants/{sample}.final.vcf"

        resistance = process_resistance(
            vcf_file, canned_text['antibiotics'], 1)

        resistance = call_resistance(
            resistance, canned_text['antibiotics'])

        result[sample] = dict()
        result[sample]['resistance'] = resistance

        all_antibiotics = (
            [antibiotic for antibiotic in resistance['resistant']] +
            [antibiotic for antibiotic in resistance['susceptible']])

        result[sample]['antibiotics'] = dict()

        for antibiotic_name in all_antibiotics:
            result[sample]['antibiotics'][antibiotic_name] = 0

        for antibiotic in resistance['resistant']:
            result[sample]['antibiotics'][antibiotic] = 1

    data = pd.DataFrame([
                {'Sample': sample, **result[sample]['antibiotics']}
                for sample in result])

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
        sample for sample in samples
        if samples[sample]['type'] == 'test_sample']

    total_test_samples = len(test_samples)

    half = math.ceil(total_test_samples/2)

    # keep count to know when to move on to next column
    count = 0

    for sample in samples:

        if samples[sample]['type'] != 'test_sample':
            continue

        qc_status = samples[sample]['qc_status']
        total_amplicons = samples[sample]['total_amplicons']
        fail_amplicons = len(samples[sample]['fail_amplicons'])
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
                    {result[sample]['resistance']['resistance_level']}
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

    for sample in samples:
        data.loc[data.samples == sample, 'text'] = (
            len(samples[sample]['passed_amplicons']))

    # TODO remove hardcoded 16 here for the total number of amplicons
    data.loc[((data.text != 16) & (data.value != 1)), 'value'] = 0.5

    mapper = LinearColorMapper(
        palette=(colors.BRAND_LIGHT_GREY, '#ffffff', colors.BRAND_BLUE),
        low=0,
        high=1
    )

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
        fill_color=transform('value', mapper)
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
    section._add_item("""
         This section summarises the antibiotic resistance assignment per
         sample. For further information on the justiciation and confidence
         of these antibiotic assignments, please review following sections
         in the report.
         <br><br>
         Interpretation:
            <ul>
            <li>
                <h5 style="display: inline;">
                    <span class="badge badge-brand-primary">
                        Resistant
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


def antibiotics_evidence(args, report_doc):
    """Retrun evidence for antibiotoci resistance."""
    json_files = args.genotype_json
    for json_file in json_files:
        with open(json_file) as json_fh:
            variant_data = json.load(json_fh)
            key = re.sub('.variants.json', '', os.path.basename(json_file))

            section = report_doc.add_section()
            section.markdown(f"###Sample: {key}")

            results = []
            for gene in variant_data:
                for antibiotic in gene['antibiotic_resistance']:
                    evidence = {}
                    evidence["gene"] = gene["gene"]
                    evidence["antibiotic"] = antibiotic["antibiotic"]
                    evidence["mutations"] = ""
                    for mutation in antibiotic['variants']:
                        if len(evidence["mutations"]) > 0:
                            evidence["mutations"] += ", "
                        evidence["mutations"] += mutation["variant"]
                        if "observed_allele" in mutation:
                            evidence["mutations"] += f"""
                             ({mutation['observed_allele']})"""
                    results.append(evidence)
            data = pd.DataFrame(results)
            placeholder = report_doc.add_section(key=f"variantSummary{key}")
            placeholder.table(
                data, index=False, key=f"variantSummary{key}",
                th_color=colors.BRAND_BLUE, paging=False, searchable=False)


def genotyped_variant_summary(samplekey, variant_data, report_doc):
    """Return genotyped variant summary."""
    results = {}
    for gene in variant_data:
        for antibiotic in gene['antibiotic_resistance']:
            for mutation in antibiotic['variants']:
                key = f'{gene["gene"]}.{mutation["variant"]}'
                if key not in results.keys():
                    triplet = mutation['triplet']
                    mutation_seq = f"{mutation['reftrip']} >" \
                        f"{triplet}" if "reftrip" in mutation else ""
                    results[key] = {
                        "gene": gene["gene"],
                        "position": mutation["position0"],
                        "variant": mutation["variant"],
                        "seqchange": mutation_seq,
                        "antibiotic": antibiotic["antibiotic"],
                        "total reads": mutation["depth"],
                        "mutant reads": mutation["support"],
                        "wild-type reads": mutation["refcount"],
                        "strand bias": mutation["fwdBias"]}
                    if ("observed_allele" in mutation and
                            mutation["observed_allele"]
                            != mutation["variant"]):
                        results[key]["variant"] += f"""
                         ({mutation['observed_allele']})
                         """
                else:
                    results[key]["antibiotic"] += f', ' \
                        f'{antibiotic["antibiotic"]}'
    data = pd.DataFrame(list(results.values()))
    if not data.empty:
        data.sort_values(by=['position'], inplace=True)
    placeholder = report_doc.add_section(key=f"variantSummary{samplekey}")
    placeholder.table(
        data, index=False, key=f"variantSummary{samplekey}",
        th_color=colors.BRAND_BLUE, paging=False, searchable=False)


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
        args, report_doc, samples, show_mapping_statistics=False):
    """Return variants evidence."""
    section = report_doc.add_section()
    section._add_item("""
        <div class="card bg-light mt-3 mb-3">
            <h5 class="card-header">
                <i class="fa fa-vials"></i> Sample Details
            </h5>
    <div class="card-body">
        The tabular information below shows a summary of
         the genetic variants observed within the sample.
         Data included in the table includes genomic coordinates
         and information on the depth of the sequence data(post-downsampling)
         and the number of reads that support the presence of the variant.
         Other supporting data includes the method used for the variant
         identification, `vcall`, and `fbias` that shows the fraction of
         reads that are found on the forward strand.
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
                    <i class="fa fa-vial"></i> {samplekey}
                </span>
                <span class="badge badge-white">
                    {samples[samplekey]['type']}
                </span>
            </h5>
            <div class="card-body">
                <h5>Variants Summary</h5>
        """)

        table = variants_table_from_vcf(vcf_file, info_fields, 1)
        placeholder = report_doc.add_section(key=f"variantSummary{samplekey}")
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
            <h5>Target Coverage by Strand</h5>
            Note that the position co-ordinate is discretized by
            {tile_size} bases.
            """)

        section.plot(target_coverage_plot_panel(
            samplekey, read_data, bed_data, ref_tiles,
            tile_size, report_doc, "Mapping data"))

        section._add_item("</div></div>")

    section._add_item("</div></div>")


def target_info(target_name, targets, reads, tile_size, ref_tiles):
    """Return info for each target."""
    print('Getting statistics for {}.'.format(target_name))

    t_reads = pr.PyRanges(reads.df[reads.df["readgroup"] == target_name])
    hits = ref_tiles.count_overlaps(t_reads, strandedness='same')

    info = targets[targets.tname == target_name].df
    # Find intersecting hits
    t_hits = hits.intersect(pr.PyRanges(info)).df
    # t_reads = rg_reads.df[rg_reads.df['readgroup'] == target_name]

    info['tsize'] = info.End - info.Start
    # use the tiles to calculate bases on target as parts of reads
    # may not overlap
    info['kbases'] = tile_size * t_hits['NumberOverlaps'].sum() / 1000
    info['median_coverage'] = t_hits.groupby(
        'Start')['NumberOverlaps'].sum().median().astype(int)
    if len(t_reads) > 0:
        info['nreads'] = t_reads.df['name'].unique().size
        info['mean_read_length'] = t_reads.df['read_length'].mean()
        info['mean_accuracy'] = t_reads.df['acc'].mean()
        fwd, rev = (len(t_reads[t_reads.df['Strand'] == x]) for x in '+-')
        info['strand_bias'] = (fwd - rev) / (fwd + rev)
    else:
        info['nreads'] = 0
        info['mean_read_length'] = np.nan
        info['mean_accuracy'] = np.nan
        info['strand_bias'] = 0
    return info


def process_samples_inputs(args):
    """Process inputs and make some decisions on QC."""
    sample_coveage = {
        sample: {
            'type': type,
            'readcount': f"{args.readcounts}/{sample}.bedtools-coverage.bed"
        } for sample, type in zip(
            args.samples, args.types
        )
    }

    controls = dict(
        test_sample=dict(threshold=args.sample_threshold),
        no_template_control=dict(threshold=args.ntc_threshold),
        positive_control=dict(threshold=args.positive_threshold))

    for sample in sample_coveage:

        sample_type = sample_coveage[sample]['type']

        coverage = process_coverage(
            sample_coveage[sample]['readcount'],
            controls[sample_type]['threshold'])

        result = determine_status(
            coverage,
            controls[sample_type]['threshold'])

        sample_coveage[sample]['total_covered'] = len(
            result['passed_targets'])

        sample_coveage[sample]['qc_status'] = result['status']
        sample_coveage[sample]['fail_amplicons'] = result['failed_targets']
        sample_coveage[sample]['passed_amplicons'] = result['passed_targets']
        sample_coveage[sample]['total_amplicons'] = (
            len(result['passed_targets'])+len(result['failed_targets']))

    return sample_coveage


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
        "--samples", nargs='+', default='unknown',
        help="git commit number")
    parser.add_argument(
        "--types", nargs='+', default='unknown',
        help="git commit number")
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

    if args.style == 'ond':
        colors = ond_colors
    elif args.style == 'ont':
        colors = ont_colors

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

    sample_types_counts = process_samples_inputs(args)

    controls(
        sample_types_counts, args.ntc_threshold,
        args.positive_threshold, report_doc)

    csv_file = 'wf-tb-amr-report.csv'
    csv_output(sample_types_counts, canned_text, csv_file)

    section_executive_summary(
        args, report_doc, sample_types_counts, canned_text)

    section_reads_per_barcode(args, report_doc)

    variants_evidence(args, report_doc, sample_types_counts)

    section = report_doc.add_section(
        section=scomponents.version_table(
            args.versions, th_color=colors.BRAND_BLUE))

    report_doc.add_section(
        section=scomponents.params_table(
            args.params, th_color=colors.BRAND_BLUE))

    report_doc.write('wf-tb-amr-report.html')


if __name__ == "__main__":
    main()
