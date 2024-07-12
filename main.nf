#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2


include { fastq_ingress } from './lib/ingress'


process getVersions {
    label "microbial"
    cpus 1
    memory "1 GB"
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


process alignReads {
    label 'microbial'
    cpus params.align_threads
    memory "12 GB"
    input:
        tuple val(sample_id), val(type), path(sample_fastq)
        path reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        tuple path("${sample_id}.bamstats"), path("${sample_id}.bam.summary"), emit: bamstats
    """
    mini_align -i ${sample_fastq} -r ${reference} -p ${sample_id}_tmp -t $task.cpus -m

    # keep only mapped reads
    samtools view --write-index -F 4 ${sample_id}_tmp.bam -o ${sample_id}.bam##idx##${sample_id}.bam.bai

    # get stats from bam
    stats_from_bam -o ${sample_id}.bamstats -s ${sample_id}.bam.summary -t $task.cpus ${sample_id}.bam
    """
}


process downSample {
    label 'microbial'
    cpus 1
    memory "1 GB"
    input:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path amplicons_bed
    output:
        tuple val(sample_id), val(type), path("${sample_id}_all_merged.sorted.bam"), path("${sample_id}_all_merged.sorted.bam.bai")
    """
    # split bam
    header_count=`samtools view -H ${sample_id}.bam | wc -l`
    lines=\$(( ${params.downsample} + \$header_count + 1 ))
    while read line;
    do
      region=`echo -e "\${line}" | cut -f1-3 | sed '1 s/\t/:/' | sed 's/\t/-/g'`;

      samtools view -bh ${sample_id}.bam \${region} > ${sample_id}_\${region}.bam;
      samtools view -h -F16 ${sample_id}_\${region}.bam > ${sample_id}_\${region}_fwd.sam;
      head -\${lines} ${sample_id}_\${region}_fwd.sam | samtools view -bh - > ${sample_id}_\${region}_fwd.bam;

      samtools view -h -f16 ${sample_id}_\${region}.bam > ${sample_id}_\${region}_rev.sam;
      head -\${lines} ${sample_id}_\${region}_rev.sam | samtools view -bh - > ${sample_id}_\${region}_rev.bam;
      samtools merge ${sample_id}_\${region}_all.bam ${sample_id}_\${region}_fwd.bam ${sample_id}_\${region}_rev.bam;

    done < ${amplicons_bed}

    samtools merge ${sample_id}_all_merged.bam *_all.bam
    samtools sort ${sample_id}_all_merged.bam > ${sample_id}_all_merged.sorted.bam
    samtools index ${sample_id}_all_merged.sorted.bam
    echo "done"

    """
}

process mpileup {
    label 'microbial'
    cpus params.mpileup_threads
    memory "12 GB"
    input:
        path reference
        path vcf_template
        path bcf_annotate_template
        tuple val(sample_id), val(type), path(bam), path(bam_index)
        path variant_db
        path variant_db_index
        path genbank
        path amplicons_bed
    output:
        tuple val(sample_id), val(type), path("${sample_id}.mpileup.annotated.processed.vcf"),  path("${sample_id}.mpileup.annotated.processed.PASS.vcf"), path(bam), path(bam_index)

    """
    #bcftools doesn't like non-normed regions
    bcftools norm -m +both ${variant_db} -Oz -o ${variant_db}.norm
    tabix ${variant_db}.norm

    # run mpileup
    bcftools mpileup \
      --max-depth 8000 \
      --threads $task.cpus \
      -BI \
      -Q 1 \
      --ff SECONDARY,UNMAP \
      --annotate INFO/AD,INFO/ADF,INFO/ADR \
      -R ${amplicons_bed} \
      -O v \
      -f ${reference} ${bam} > ${sample_id}.mpileup.vcf

    bgzip ${sample_id}.mpileup.vcf
    tabix ${sample_id}.mpileup.vcf.gz

    bcftools norm --remove-duplicates -Oz ${sample_id}.mpileup.vcf.gz -o ${sample_id}.mpileup.vcf.gz.dedup
    tabix ${sample_id}.mpileup.vcf.gz.dedup

    bcftools norm -m- -Oz ${sample_id}.mpileup.vcf.gz.dedup -o ${sample_id}.mpileup.vcf.gz.norm
    tabix ${sample_id}.mpileup.vcf.gz.norm

    bcftools view ${sample_id}.mpileup.vcf.gz.norm > ${sample_id}.mpileup.vcf.gz.norm.vcf


    # call variants from pileup
    workflow-glue process_mpileup \
      --template ${vcf_template} \
      --mpileup ${sample_id}.mpileup.vcf.gz.norm.vcf \
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
    cpus 4
    memory "12 GB"
    input:
        path reference
        path genbank
        path variant_db
        path vcf_template
        path bcf_annotate_template
        tuple val(sample_id), val(type), path("${sample_id}.mpileup.annotated.processed.vcf"),  path("${sample_id}.mpileup.annotated.processed.PASS.vcf"), path(bam), path(bam_index)
    output:
        tuple val(sample_id), val(type), path("${sample_id}.final.vcf")

    """
    # whatshap needs read group
    samtools addreplacerg -r "ID:${sample_id}\tSM:${sample_id}" -o ${sample_id}.rg.bam ${bam}
    samtools index ${sample_id}.rg.bam

    # index fasta
    samtools faidx ${reference}

    read_count=`samtools view -c ${sample_id}.rg.bam`

    if [ "\${read_count}" -gt "0" ]; then

      # phase variants
      whatshap phase \
        -o ${sample_id}.phased.vcf \
        --reference=${reference} \
        ${sample_id}.mpileup.annotated.processed.PASS.vcf ${sample_id}.rg.bam
    else
      cp ${vcf_template} ${sample_id}.phased.vcf
    fi

    # add codon numbers to those variants which are not in our db but we want to phase because they could affect the same codon
    vcf-annotator ${sample_id}.phased.vcf ${genbank} > ${sample_id}.phased.codon.vcf

    # process phased variants
    workflow-glue process_whatshap \
      --phased_vcf ${sample_id}.phased.codon.vcf \
      --out_vcf ${sample_id}.phased.processed.vcf \
      --template ${vcf_template} \
      --sample ${sample_id}

    # sort
    bcftools sort ${sample_id}.phased.processed.vcf > ${sample_id}.phased.processed.sorted.vcf

    # re-annotate our newly phased variants
    bgzip ${sample_id}.phased.processed.sorted.vcf
    tabix ${sample_id}.phased.processed.sorted.vcf.gz
    tabix ${variant_db}

    bcftools annotate \
      -c CHROM,POS,REF,GENE,STRAND,AA,FEATURE_TYPE,EFFECT,GENE_LOCUS,WHO_POS,ANTIBIOTICS,PROTEIN_ID,HGVS_NUCLEOTIDE,HGVS_PROTEIN,CODON_NUMBER,ORIGIN \
      --remove INFO/FeatureType,INFO/IsSynonymous,INFO/IsTransition,INFO/IsGenic,INFO/IsPseudo,INFO/Inference,INFO/AltCodon,INFO/AltAminoAcid,INFO/Note,INFO/AminoAcidChange,INFO/Product,INFO/SNPCodonPosition \
      -h ${bcf_annotate_template} \
      -a ${variant_db} \
      ${sample_id}.phased.processed.sorted.vcf.gz > ${sample_id}.phased.processed.sorted.annotated.vcf

    # filter out those without annotation - they are not in the WHO database
    bcftools filter -i 'INFO/ORIGIN=="WHO_CANONICAL"' ${sample_id}.phased.processed.sorted.annotated.vcf > ${sample_id}.final.vcf
    """
}


process countReadsRegions {
    label "microbial"
    cpus 1
    memory "1 GB"
    input:
        path amplicons_bed
        tuple val(sample_id), val(type), path(bam), path(bai)
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bedtools-coverage.bed"), emit: bed_files

    """
    samtools view -q 1 -bh ${bam} | bedtools coverage -d -a ${amplicons_bed} -b - > ${sample_id}.bedtools-coverage.bed
    """
}


process report {
    label "microbial"
    cpus 1
    memory "2 GB"
    input:
        val metadata
        path "bed_files/*"
        path "per_barcode_stats/?.gz"
        path "variants/*"
        path "params.json"
        path "pickedreads/*"
        path reference
        path "versions/*"
        path amplicons_bed
        path report_config
    output:
        tuple path("wf-tb-amr-report.html"), path("wf-tb-amr-report.csv"), emit: run_report
        path("jsons/*"), emit: jsons
    script:
    def metadata = new JsonBuilder(metadata).toPrettyString()

    myDir = file("jsons")
    myDir.mkdirs()
    """
    echo '${metadata}' > metadata.json
    mkdir jsons
    workflow-glue report \
        --revision $workflow.revision \
        --commit $workflow.commitId \
        --per_barcode_stats per_barcode_stats/* \
        --bed ${amplicons_bed} \
        --genotype_json variants/* \
        --params params.json \
        --pickedreads pickedreads/* \
        --reference $reference \
        --metadata metadata.json \
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
    memory "2 GB"
    input:
      val samples
      path "jsons/*"
      path report_config
    output:
      path "${samples.sample_id[0]}_report.html"
    """
    workflow-glue report_single_sample \
      --revision $workflow.revision \
      --commit $workflow.commitId \
      --canned_text ${report_config} \
      --sample_id ${samples.sample_id[0]} \
      --barcode ${samples.barcode[0]} \
      --jsons jsons/* \
      --ntc_threshold="${params.ntc_threshold}" \
      --pos_threshold="${params.positive_threshold}" \
      --sample_threshold="${params.sample_threshold}" \
      --group 1 \
      --style ont
    """
}

process getParams {
    label "microbial"
    cpus 1
    memory "1 GB"
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
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
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
        samples_for_processing = samples.map {it -> [it[0].alias, it[0].type, it[1]]}
        fastcat_stats = samples.map {it -> it[2]}
        
        software_versions = getVersions()

        workflow_params = getParams()

        // do alignment
        alignments = alignReads(samples_for_processing, reference)

        // do crude downsampling
        if (params.downsample != null){
          log.warn("Downsampling data to ${params.downsample} in each direction.")
          downsample = downSample(alignments[0], amplicons_bed)
        } else {
          log.warn("Using all data, downsampling level has not been specified.")
          downsample = alignments
        }

        // do mpileup
        mpileup_result = mpileup(reference, vcf_template, bcf_annotate_template, downsample[0], variant_db, variant_db+".tbi", genbank, amplicons_bed)

        // phase variants
        whatshap_result = whatshap(reference, genbank, variant_db, vcf_template, bcf_annotate_template, mpileup_result)

        // do some coverage calcs
        region_read_count = countReadsRegions(amplicons_bed, downsample[0])

        // we need a channel from the sample sheet in case the NTC has no reads
        sample_sheet = Channel.fromPath(params.sample_sheet).splitCsv(skip: 1, sep: ",").map{ it -> tuple("sample_id":it[1], "type":it[2], "barcode":it[0]) }

        // generate run report
      	report = report(
              sample_sheet.collect(),
              region_read_count.bed_files.map{it[2]}.collect(),
              samples.map { it[2].resolve("per-read-stats.tsv.gz") }.toList(),
              whatshap_result.map{it[2]}.collect(),
      	      workflow_params,
              alignments.bamstats.collect(),
              reference,
              software_versions.collect(),
              amplicons_bed,
              report_config
      	)

        test_samples = sample_sheet.filter{it[0].type == "test_sample"}

        // Generate single sample report
        report_single_sample = reportSingle(
              test_samples,
              report.jsons.collect(),
              report_config
        )

        output_alignments = alignments[0].map{ it -> return tuple(it[2], it[3]) }

        results = report.run_report.mix(
            output_alignments.collect(),
            whatshap_result.map{ it[2]}.collect(),
            report_single_sample.collect()
        )


    emit:
        results
        telemetry = workflow_params
}



// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

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

    //filter unclassified here
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "stats":true,
        "per_read_stats":true ])
        // "sample_sheet":params.sample_sheet]).filter { it[1].sample_id != "unclassified" }

      //get reference
    if (params.reference == null){
      params.remove('reference')
      params._reference = projectDir.resolve("./data/primer_schemes/V3/NC_000962.3.fasta").toString()
    } else {
      params._reference = file(params.reference, type: "file", checkIfExists:true).toString()
      params.remove('reference')
    }

    // Variant DB
    if (params.variant_db == null){
      params.remove('variant_db')
      params._variant_db = projectDir.resolve("./data/primer_schemes/V3/variant_db.sorted.normalised.vcf.gz").toString()
    } else {
      params._variant_db = file(params.variant_db, type: "file", checkIfExists:true).toString()
      params.remove('variant_db')
    }

    // Genbank
    if (params.genbank == null){
      params.remove('genbank')
      params._genbank = projectDir.resolve("./data/primer_schemes/V3/NC_000962.3.gb").toString()
    } else {
      params._genbank = file(params.genbank, type: "file", checkIfExists:true).toString()
      params.remove('genbank')
    }

    // TB amplicons
    if (params.amplicons_bed == null){
      params.remove('amplicons_bed')
      params._amplicons_bed = projectDir.resolve("./data/primer_schemes/V3/TB_amplicons.bed").toString()
    } else {
      params._amplicons_bed = file(params.reference, type: "file", checkIfExists:true).toString()
      params.remove('amplicons_bed')
    }

    // Single sample report text
    if (params.report_config == null){
      params.remove('report_config')
      params._report_config = projectDir.resolve("./data/primer_schemes/V3/report_config.eng.json").toString()
    } else {
      params._report_config = file(params.report_config, type: "file", checkIfExists:true).toString()
      params.remove('report_config')
    }

    vcf_template = projectDir.resolve("./data/template.vcf").toString()
    bcf_annotate_template = projectDir.resolve("./data/bcftools_annotate_header.txt").toString()

    pipeline(samples, file(params._reference), file(params._amplicons_bed), file(params._variant_db), file(params._genbank), file(vcf_template), file(bcf_annotate_template), file(params._report_config))
    
    output(pipeline.out.results)
}


workflow.onComplete {
  Pinguscript.ping_complete(nextflow, workflow, params)
}

workflow.onError {
  Pinguscript.ping_error(nextflow, workflow, params)
}

