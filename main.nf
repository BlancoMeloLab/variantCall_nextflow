#!/usr/bin/env nextflow
/*
========================================================================================
                         SteadyFlow - Variant Calling PIPELINE
========================================================================================
Steady Sate Analysis Pipeline. Started 2022-10-17.
 #### Homepage / Documentation

 # Author: Qing Yang <qing.yang@colorado.edu>

========================================================================================
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq
        1a. SRA tools -- fasterq-dump sra to generate fastq file
        1b. FastQC (pre-trim) -- perform pre-trim FastQC on fastq files

    2. Trimming
        2a. BBDuk -- trim fastq files for quality and adapters
        2b. FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)

    3. Mapping w/ BWA -- map to genome reference file and merge bam files

    4. Mark duplicates and base quality calibration

    5. Haplotype calling


*/


def helpMessage() {
    log.info"""
    =========================================
     SteadyFlow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile slurm --fastqs '/project/*_{R1,R2}*.fastq' --outdir '/project/'

    Required arguments:
         -profile                      Configuration profile to use. <base, slurm>
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         --email                       Where to send workflow report email.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --flip                         Reverse complements each strand. Necessary for some library preps.
        --flipR2                       Reverse complements R2 only.

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.

    QC Options:
        --skipMultiQC                  Skip running MultiQC.

    Analysis Options:
        --MuTect2
        --HaplotypeCaller

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.bbmap_adapters = "$baseDir/bin/adapters.fa"
params.workdir = "./nextflowTemp"

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.chrom_sizes ){
    chrom_sizes = file(params.chrom_sizes)
    if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }

 if ( params.picard_path ){
     picard_path = file("${params.picard_path}")
 }
 if ( params.known_sites_dir ){
     known_sites_dir = file("${params.known_sites_dir}")
 }

if ( params.bwa_indices ){
    bwa_indices = file("${params.bwa_indices}")
}

if ( params.bbmap_adapters){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * Create a channel for input read files
 */
if (params.fastqs) {
    if (params.singleEnd) {
        fastq_reads_qc = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
    } else {
        Channel
            .fromFilePairs( params.fastqs, size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
    }
}
else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
}

if (params.sras) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
}
else {
    read_files_sra = Channel.empty()
}


// Header log info
log.info """=======================================================
RNAseq Flow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'variantCallFlow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.reads) summary['Reads']     = params.reads
if(params.fastqs) summary['Fastqs']   = params.fastqs
if(params.sras) summary['SRAs']       = params.sras
summary['Genome Ref']       = params.genome
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Reverse Comp']     = params.flip ? 'YES' : 'NO'
summary['Reverse Comp R2']  = params.flipR2 ? 'YES' : 'NO'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/software_versions/", mode: 'copy', pattern: '*.txt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    multiqc --version > v_multiqc.txt

    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1 -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"
    cpus 1
    if (params.savefq || params.saveAllfq) {
        publishDir "${params.outdir}/fastq", mode: 'copy'
    }

    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("*.fastq.gz") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_gzip_sra


    script:
    prefix = reads.baseName
    if (!params.singleEnd) {
        """
        echo ${prefix}

        fastq-dump --split-3 ${reads} --gzip
        """
    } else {
         """
        fastq-dump ${reads} --gzip
        """
    }

}

/*
 * STEP 2 - FastQC
 */

process fastQC {
    errorStrategy { task.exitStatus=0 ? 'ignore' : 'terminate' }

    tag "$prefix"
    memory '8 GB'
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(prefix), file(reads) from fastq_reads_qc.mix(fastq_reads_qc_sra)

    output:
    file "*.{zip,html,txt}" into fastqc_results

    script:
    """
    fastqc $reads
    """
}


/*
 * STEP 3 - Trimming
 */

process bbduk {
    errorStrategy { task.exitStatus=0 ? 'ignore' : 'terminate' }

    tag "$name"
    cpus 16
    time '2h'
    memory '20 GB'
    publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*.txt"
    if (params.saveTrim || params.saveAllfq) {
        publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    }

    input:
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    set val(name), file ("*.trim.fastq.gz") into trimmed_reads_fastqc, trimmed_reads_bwa
    file "*.txt" into trim_stats

    script:
    if (!params.singleEnd && params.flip) {
        """
        echo ${name}

        reformat.sh -Xmx20g \
                t=16 \
                in=${name}_1.fastq.gz \
                in2=${name}_2.fastq.gz \
                out=${name}_1.flip.fastq.gz \
                out2=${name}_2.flip.fastq.gz \
                rcomp=t

        bbduk.sh -Xmx20g \
                t=16 \
                in=${name}_1.flip.fastq.gz \
                in2=${name}_2.flip.fastq.gz \
                out=${name}_R1.flip.trim.fastq.gz \
                out2=${name}_R2.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=21 mink=11 hdist=2 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                stats=${name}.trimstats.txt \
                refstats=${name}.refstats.txt
        """
    } else if (!params.singleEnd && params.flipR2) {
                """
        echo ${name}

        reformat.sh -Xmx20g \
                t=16 \
                in=${name}_1.fastq.gz \
                in2=${name}_2.fastq.gz \
                out=${name}_1.flip.fastq.gz \
                out2=${name}_2.flip.fastq.gz \
                rcompmate=t

        bbduk.sh -Xmx20g \
                t=16 \
                in=${name}_1.flip.fastq.gz \
                in2=${name}_2.flip.fastq.gz \
                out=${name}_R1.flip.trim.fastq.gz \
                out2=${name}_R2.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=21 mink=11 hdist=2 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                stats=${name}.trimstats.txt \
                refstats=${name}.refstats.txt
        """
    } else if (params.flip) {
        """
        echo ${name}


        reformat.sh -Xmx20g \
                t=16 \
                in=${name}.fastq.gz \
                out=${name}.flip.fastq.gz \
                rcomp=t


        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}.flip.fastq.gz \
                  out=${name}.flip.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=21 mink=11 hdist=2 \
                  nullifybrokenquality=t \
                  maq=10 minlen=25 \
                  tpe tbo \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    }
        else if (!params.singleEnd) {
        """
        echo ${name}

        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}_1.fastq.gz \
                  in2=${name}_2.fastq.gz \
                  out=${name}_R1.trim.fastq.gz \
                  out2=${name}_R2.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    } else {
        """
        echo ${name}

        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}.fastq.gz \
                  out=${name}.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    }
}


/*
 * STEP 4 - Trimmed FastQC
 */

process fastqc_trimmed {
    tag "$name"
    memory '4 GB'
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html,txt}" into trimmed_fastqc_results

    script:
    """

    echo ${name}

    fastqc $trimmed_reads
    """
}


/*
 * STEP 5 - Map reads to reference genome
 */

process bwa {
    errorStrategy { task.exitStatus=0 ? 'ignore' : 'terminate' }
    tag "$name"
    cpus 16
    memory '50 GB'
    time '5h'

    input:
    val(indices_path) from bwa_indices
    val(picard_path) from picard_path
    set val(name), file(trimmed_reads) from trimmed_reads_bwa

    output:
    set val(name), file("*.bam") into bwa_bam

    script:
    if (!params.singleEnd) {
        """
        echo ${name}
        bwa mem -t 16 \
                -T 0 \
                ${indices_path} \
                ${name}_R1.trim.fastq.gz ${name}_R2.trim.fastq.gz | \
                java -Xmx30g -jar ${picard_path}/picard.jar \
                SortSam \
                I=/dev/stdin \
                O=${name}.bam \
                SORT_ORDER=coordinate

        """
    } else {
        """
        echo ${name}
        bwa mem -t 16 \
                -T 0 \
                ${indices_path} \
                ${trimmed_reads} | \
                java -Xmx30g -jar ${picard_path}/picard.jar \
                SortSam \
                I=/dev/stdin \
                O=${name}.bam \
                SORT_ORDER=coordinate

        """
    }
}

/*
 * STEP 6 - Mark duplicates, base quality calibration
 */

process preprocess {
    tag "$name"
    memory '100 GB'
    cpus 16
    publishDir "${params.outdir}" , mode: 'copy',
                saveAs: {filename ->
                  if (filename.indexOf("preprocessed.bam") > 0) "mapped/bams/$filename"
                }

    input:
    val(known_sites_dir) from known_sites_dir
    val(genome) from genome
    val(picard_path) from picard_path
    set val(name), file(mapped_bam) from bwa_bam

    output:
    set val(name), file("${name}.preprocessed.bam") into preproc_bam
    file("metrics.txt") into bam_metrics

    script:
    """
    gatk MarkDuplicates \
          --java-options -Xmx30g \
          I=${name}.bam \
          O=${name}.marke_dup.bam \
          M=metrics.txt

    java -Xmx30g -jar ${picard_path}/picard.jar AddOrReplaceReadGroups \
           I=${name}.marke_dup.bam \
           O=${name}.rg.bam \
           RGID=${name} \
           RGLB=lib1 \
           RGPL=ILLUMINA \
           RGPU=unit1 \
           RGSM=${name}

    gatk BaseRecalibrator \
          --java-options -Xmx30g \
          --input ${name}.rg.bam \
          --output ${name}.recal.table \
          --known-sites ${known_sites_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
          --known-sites ${known_sites_dir}/Homo_sapiens_assembly38.dbsnp138.vcf \
          --known-sites ${known_sites_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
          --reference ${genome}

    gatk ApplyBQSR \
      -R ${genome} \
      -I ${name}.rg.bam \
      --bqsr-recal-file ${name}.recal.table \
      -O ${name}.preprocessed.bam

    """
}

/*
 * WIP STEP 7 - Individual variant calls, run in GVCF format
 */
 process haplotypeCaller {
     tag "$name"
     memory '100 GB'
     cpus 16
     publishDir "${params.outdir}/processed/vcf/", mode: 'copy',
         saveAs: {filename -> if(filename.indexOf("vcf") > 0) "$filename"}

     input:
     val(genome) from genome
     set val(name), file(preproc_bam) from preproc_bam

     output:
     set val(name), file("${name}.g.vcf.gz") into gvcf

     script:
     """
     gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R ${genome} \
      -I ${preproc_bam} \
      -O ${name}.g.vcf.gz \
      -ERC GVCF
     """

 }

/*
 * STEP 9 - MultiQC
 */
process multiQC {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC

    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}


/*
 * STEP 10 - Output Description HTML
 */
//
//process output_documentation {
//    tag "$prefix"
//    publishDir "${params.outdir}/pipeline_info/", mode: 'copy'
//
//    input:
//    file output_docs
//
//    output:
//    file "results_description.html"
//
//    script:
//    """
//    markdown_to_html.r $output_docs results_description.html
//    """
//}
//


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[SteadyFlow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[SteadyFlow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[SteadyFlow] Pipeline Complete"

}
