#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2


include { fastq_ingress } from './lib/fastqingress'
include { start_ping; end_ping } from './lib/ping'


process getVersions {
    label "microbial"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    bcftools --version | grep bcftools | sed 's/ /,/' >> versions.txt
    """
}


process combineFastq {
    label 'microbial'
    cpus 1
    input:
        tuple val(sample_id), val(barcode), file(directory), val(type)
    output:
        tuple val(sample_id), val(type), path("${sample_id}.fastq.gz"), emit: sample
        path "${sample_id}.stats", emit: fastqstats
    """
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} | seqkit seq -m 200 - > ${sample_id}.fastq
    gzip ${sample_id}.fastq
    """
}


process alignReads {
    label 'microbial'
    cpus params.threads
    input:
        tuple val(sample_id), val(type), file(sample_fastq)
        file reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        tuple path("${sample_id}.bamstats"), path("${sample_id}.bam.summary"), emit: bamstats
    """
    mini_align -i ${sample_fastq} -r ${reference} -p ${sample_id} -t $task.cpus -m
    stats_from_bam -o ${sample_id}.bamstats -s ${sample_id}.bam.summary -t $task.cpus ${sample_id}.bam
    """
}


process mpileup {
    label 'microbial'
    cpus params.threads
    input:
        file reference
        file vcf_template
        file bcf_annotate_template
        tuple val(sample_id), val(type), path(bam), path(bam_index)
        file variant_db
        file varinat_db_index
        file genbank
    output:
        tuple val(sample_id), val(type), path("${sample_id}.mpileup.annotated.processed.vcf"),  path("${sample_id}.mpileup.annotated.processed.PASS.vcf"), path(bam), path(bam_index)

    """
    #bcftools doesn't like non-normed regions
    bcftools norm -m +both ${variant_db} -Oz -o ${variant_db}.norm
    tabix ${variant_db}.norm

    # run mpileup
    bcftools mpileup \
      --max-depth 8000 \
      --threads ${params.threads} \
      -B \
      -Q 1 \
      --ff SECONDARY,UNMAP \
      --annotate INFO/AD,INFO/ADF,INFO/ADR \
      -R ${variant_db}.norm \
      -O v \
      -f ${reference} ${bam} > ${sample_id}.mpileup.vcf

    bgzip ${sample_id}.mpileup.vcf
    tabix ${sample_id}.mpileup.vcf.gz

    bcftools norm --remove-duplicates -Oz ${sample_id}.mpileup.vcf.gz -o ${sample_id}.mpileup.vcf.gz.dedup
    tabix ${sample_id}.mpileup.vcf.gz.dedup

    bcftools norm -m- -Oz ${sample_id}.mpileup.vcf.gz.dedup -o ${sample_id}.mpileup.vcf.gz.norm
    tabix ${sample_id}.mpileup.vcf.gz.norm

    # annotate pileup - if you don't use non-norm db then you lose annotations
    bcftools annotate \
      -c CHROM,POS,REF,GENE,STRAND,AA,FEATURE_TYPE,EFFECT,GENE_LOCUS,WHO_POS,ANTIBIOTICS,PROTEIN_ID,HGVS_NUCLEOTIDE,HGVS_PROTEIN,CODON_NUMBER,ORIGIN \
      -h ${bcf_annotate_template} \
      -a ${variant_db} \
      ${sample_id}.mpileup.vcf.gz.norm > ${sample_id}.mpileup.annotated.vcf

    # call variants from pileup
    process_mpileup.py \
      --template ${vcf_template} \
      --mpileup ${sample_id}.mpileup.annotated.vcf \
      --out_vcf ${sample_id}.mpileup.annotated.processed.vcf \
      --sample ${sample_id} \
      -a $params.maf \
      -d $params.minimum_read_support \
      -b $params.strand_bias \
      -p 20

    # filter PASS variants
    bcftools view --exclude-type indels ${sample_id}.mpileup.annotated.processed.vcf | bcftools view -f 'PASS' - > ${sample_id}.mpileup.annotated.processed.PASS.vcf
    """

}


process whatshap {
    label 'microbial'
    cpus 1
    input:
        file reference
        file genbank
        file variant_db
        file vcf_template
        file bcf_annotate_template
        tuple val(sample_id), val(type), path("${sample_id}.mpileup.annotated.processed.vcf"),  path("${sample_id}.mpileup.annotated.processed.PASS.vcf"), path(bam), path(bam_index)
    output:
        tuple val(sample_id), val(type), path("${sample_id}.final.vcf")

    """
    # whatshap needs read group
    samtools addreplacerg -r "ID:${sample_id}\tSM:${sample_id}" -o ${sample_id}.rg.bam ${bam}
    samtools index ${sample_id}.rg.bam

    # index fasta
    samtools faidx ${reference}

    # phase variants
    whatshap phase \
      -o ${sample_id}.phased.vcf \
      --reference=${reference} \
      ${sample_id}.mpileup.annotated.processed.PASS.vcf ${sample_id}.rg.bam

    # add codon numbers to those variants which are noit in our db but we want to phase because they could affect the same codon
    vcf-annotator ${sample_id}.phased.vcf ${genbank} > ${sample_id}.phased.codon.vcf

    # process phased variants
    process_whatshap.py \
      --phased_vcf ${sample_id}.phased.codon.vcf \
      --out_vcf ${sample_id}.phased.processed.vcf \
      --template ${vcf_template}

    # sort
    bcftools sort ${sample_id}.phased.processed.vcf > ${sample_id}.phased.processed.sorted.vcf

    # re-annotate our newly phased variants
    bgzip ${sample_id}.phased.processed.sorted.vcf
    tabix ${sample_id}.phased.processed.sorted.vcf.gz
    tabix ${variant_db}

    bcftools annotate \
      -c CHROM,POS,REF,GENE,STRAND,AA,FEATURE_TYPE,EFFECT,GENE_LOCUS,WHO_POS,ANTIBIOTICS,PROTEIN_ID,HGVS_NUCLEOTIDE,HGVS_PROTEIN,CODON_NUMBER,ORIGIN \
      -h ${bcf_annotate_template} \
      -a ${variant_db} \
      ${sample_id}.phased.processed.sorted.vcf.gz > ${sample_id}.phased.processed.sorted.annotated.vcf

    # filter out those without annotation - they are not in the WHO database
    bcftools filter -i 'INFO/ORIGIN=="WHO_CANONICAL"' ${sample_id}.phased.processed.sorted.annotated.vcf > ${sample_id}.final.vcf
    """
}


process countReadsRegions {
    label "microbial"
    input:
        tuple val(sample_id), val(type), path(bam), path(bai)
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bedtools-coverage.bed"), emit: bed_files

    """
    samtools view -q 1 -bh ${bam} | bedtools coverage -d -a ${params._amplicons_bed} -b - > ${sample_id}.bedtools-coverage.bed
    """
}


process report {
    label "microbial"
    cpus 1
    input:
        val samples
        val types
        path "bed_files/*"
	      path "per_barcode_stats/*"
        path "variants/*"
        path "params.json"
        path "pickedreads/*"
        file reference
        path "versions/*"
        path amplicons_bed
        file report_config
    output:
        path "*report.html", emit: html
    """
    report.py \
        --revision $workflow.revision \
        --commit $workflow.commitId \
        --per_barcode_stats per_barcode_stats/* \
        --bed ${amplicons_bed} \
        --genotype_json variants/* \
        --params params.json \
        --pickedreads pickedreads/* \
        --reference $reference \
        --samples $samples \
        --types $types \
        --readcounts bed_files \
        --ntc_threshold="${params.ntc_threshold}" \
        --positive_threshold="${params.positive_threshold}" \
        --sample_threshold="${params.sample_threshold}" \
        --versions versions \
        --canned_text ${report_config} \
        --style ont
    """
}


process reportSingle {
    label "microbial"
    cpus 1
    input:
      tuple val(sample_id), val(barcode), path(fastq), val(sample_type), val(sample_type), file(variants), val(sample_type), file(coverage)
      file ntc_coverage
      file pos_coverage
      file report_config
    output:
      file "${sample_id}_report.html"
    """
    report_single_sample.py \
      --revision $workflow.revision \
      --commit $workflow.commitId \
      --canned_text ${report_config} \
      --vcf ${variants} \
      --sample_id ${sample_id} \
      --barcode ${barcode} \
      --pos_coverage_file ${pos_coverage} \
      --ntc_coverage_file ${ntc_coverage} \
      --sample_coverage_file ${coverage} \
      --ntc_threshold="${params.ntc_threshold}" \
      --pos_threshold="${params.positive_threshold}" \
      --sample_threshold="${params.sample_threshold}" \
      --group 1 \
      --style ont
    """
}


process reportAppendix {
    label "microbial"
    cpus 1
    input:
      tuple val(sample_id), val(barcode), path(fastq), val(sample_type), val(sample_type), file(variants), val(sample_type), file(coverage)
      file ntc_coverage
      file pos_coverage
      file report_config
    output:
      file "${sample_id}_appendix.html"
    """
    report_appendix.py \
      --revision $workflow.revision \
      --commit $workflow.commitId \
      --canned_text ${report_config} \
      --vcf ${variants} \
      --sample_id ${sample_id} \
      --barcode ${barcode} \
      --pos_coverage_file ${pos_coverage} \
      --ntc_coverage_file ${ntc_coverage} \
      --sample_coverage_file ${coverage} \
      --ntc_threshold="${params.ntc_threshold}" \
      --pos_threshold="${params.positive_threshold}" \
      --sample_threshold="${params.sample_threshold}" \
      --group 2 3 \
      --style ont
    """
}


process getParams {
    label "microbial"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "microbial"

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """

}


workflow pipeline {
    take:
        samples
        reference
        amplicons_bed
        variant_db
        genbank
        vcf_template
        bcf_annotate_template
        report_config
    main:

        software_versions = getVersions()

        workflow_params = getParams()

        // combine fastq files
        sample_fastqs = combineFastq(samples)

        // do alignment
        alignments = alignReads(sample_fastqs.sample, reference)

        // do mpileup
        mpileup_result = mpileup(reference, vcf_template, bcf_annotate_template, alignments[0], variant_db, variant_db+".tbi", genbank)

        // phase variants
        whatshap_result = whatshap(reference, genbank, variant_db, vcf_template, bcf_annotate_template, mpileup_result)

        // do some coverage calcs
        region_read_count = countReadsRegions(alignments[0])

        samples_region = region_read_count.bed_files.map{it[0]}.collect().map{it.join(' ')}
        types = region_read_count.bed_files.map{it[1]}.collect().map{it.join(' ')}

        // generate run report
      	report = report(
              samples_region,
              types,
              region_read_count.bed_files.map{it[2]}.collect(),
      	      sample_fastqs.fastqstats.collect(),
              whatshap_result.map{it[2]}.collect(),
      	      workflow_params,
              alignments.bamstats.collect(),
              reference,
              software_versions.collect(),
              amplicons_bed,
              report_config
      	)

        // get barcodes for the called_variants channel
        for_report = samples.join(whatshap_result.join(region_read_count))

        // we want to deal with each sample with the controls so that
        // we can decide on the validity of the result
        test_samples = for_report.filter {it[3] == "test_sample"}
        ntc_samples = for_report.filter {it[3] == "no_template_control"}
        positive_samples = for_report.filter {it[3] == "positive_control"}

        // Generate single sample report
        report_single_sample = reportSingle(
              test_samples,
              ntc_samples.map{ it[7] }.collect(),
              positive_samples.map{ it[7] }.collect(),
              report_config
        )

        // Generate additional single sample report on wider set of variants
        report_appendix = reportAppendix(
              test_samples,
              ntc_samples.map{ it[7] }.collect(),
              positive_samples.map{ it[7] }.collect(),
              report_config
        )

        output_alignments = alignments[0].map{ it -> return tuple(it[2], it[3]) }


        samples = region_read_count.map{ it[0]}.collect().map{ it.join(' ')}
        types = region_read_count.map{ it[1]}.collect().map{ it.join(' ')}
        bed_files = region_read_count.map{ it[2]}.collect().map{ it.join(' ')}


        results = report.concat(
            whatshap_result.map{ it[2]}.collect(),
            whatshap_result.collect(),
            output_alignments.collect(),
            report_single_sample.collect(),
            report_appendix.collect()
        )


    emit:
        results
        telemetry = workflow_params
}



// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }

    start_ping()
    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)

      //get reference
      if (params.reference == null){
        params.remove('reference')
        params._reference = projectDir.resolve("./data/primer_schemes/V2/NC_000962.3.fasta").toString()
      } else {
        params._reference = file(params.reference, type: "file", checkIfExists:true).toString()
        params.remove('reference')
      }

      // Variant DB
      if (params.variant_db == null){
        params.remove('variant_db')
        params._variant_db = projectDir.resolve("./data/primer_schemes/V2/variant_db.sorted.normalised.vcf.gz").toString()
      } else {
        params._variant_db = file(params.variant_db, type: "file", checkIfExists:true).toString()
        params.remove('variant_db')
      }

      // Genbank
      if (params.genbank == null){
        params.remove('genbank')
        params._genbank = projectDir.resolve("./data/primer_schemes/V2/NC_000962.3.gb").toString()
      } else {
        params._genbank = file(params.genbank, type: "file", checkIfExists:true).toString()
        params.remove('genbank')
      }

      // TB amplicons
      if (params.amplicons_bed == null){
        params.remove('amplicons_bed')
        params._amplicons_bed = projectDir.resolve("./data/primer_schemes/V2/TB_amplicons.bed").toString()
      } else {
        params._amplicons_bed = file(params.reference, type: "file", checkIfExists:true).toString()
        params.remove('amplicons_bed')
      }

      // Single sample report text
      if (params.report_config == null){
        params.remove('report_config')
        params._report_config = projectDir.resolve("./data/primer_schemes/V2/report_config.eng.json").toString()
      } else {
        params._report_config = file(params.report_config, type: "file", checkIfExists:true).toString()
        params.remove('report_config')
      }

    vcf_template = projectDir.resolve("./data/template.vcf").toString()
    bcf_annotate_template = projectDir.resolve("./data/bcftools_annotate_header.txt").toString()

    pipeline(samples, file(params._reference), file(params._amplicons_bed), file(params._variant_db), file(params._genbank), file(vcf_template), file(bcf_annotate_template), file(params._report_config))

    output(pipeline.out.results)

    end_ping(pipeline.out.telemetry)
}
