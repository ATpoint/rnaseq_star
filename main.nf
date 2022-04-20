#! /usr/bin/env nextflow

/*
   Starting from paired-end RNA-seq data: 
   Run fastqc, align to existing index with STAR,
   return unsorted BAM files, run featureCounts to get a count matrix.
*/

nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Intro
//--------------------------------------------------------------------------

ANSI_RESET   = "\u001B[0m"
ANSI_RED     = "\u001B[31m"
ANSI_YELLOW  = "\u001B[33m"
ANSI_GREEN   = "\u001B[32m"
DASHEDDOUBLE = "=".multiply(120)

Date date = new Date()
String datePart = date.format("yyyy-dd-MM -- ")
String timePart = date.format("HH:mm:ss")
def start_date = datePart + timePart

println "$ANSI_GREEN" + "$DASHEDDOUBLE"
println "Pipeline:      rnaseq_star"
println "GitHub:        https://github.com/ATpoint/rnaseq_star/"
println "Author:        Alexander Toenges (@ATpoint)"
println "Runname:       $workflow.runName"
println "Profile:       $workflow.profile"
println "Start:         $start_date"
println ""
println "Params Summary:"

def max_char = params.keySet().collect { it.length() }.max()  
params.each { name, entry -> 

    // eye-candy, what (or not) to include into the params summary
    if(name=="genome" | 
       name=="run_index" | 
       name=="is_slurm" | 
       name=="is_container" | 
       !entry.toString()?.trim()) return

    if(!params.is_slurm & name=="clusteropts") return
    if(!params.is_slurm & name=="queue") return
    if(!params.is_container & name=="container") return

    def use_length = max_char - name.length()
    def spacer = ' '.multiply(use_length)
    println "${name} ${spacer}:: ${entry}" 

}
println "${DASHEDDOUBLE}"
println "$ANSI_RESET"

//--------------------------------------------------------------------------
// Processes
//--------------------------------------------------------------------------

process FASTQC {

    tag "$sample_id"

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir "${params.outdir}/fastqc/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
            
    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    
    script: 
    """
    fastqc --threads 1 $reads[0]
    fastqc --threads 1 $reads[1]
    """     

}

// this is only for the CI tests as STAR index is too large to fit into GitHub as an asset
process INDEX {

    cpus 1
    memory 1.GB

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(genome)
    path(gtf)

    output:
    path "index", emit: idx

    script:
    """
    mkdir index
    STAR \
        --runMode genomeGenerate \
        --runThreadN $task.cpus \
        --sjdbOverhang 149 \
        --genomeDir index \
        --genomeFastaFiles $genome \
        --sjdbGTFfile $gtf
    """

}

process STAR {

    tag "$sample_id"
    
    errorStrategy 'finish'

    publishDir "${params.outdir}/star/", mode: 'copy'
    
    cpus params.cpus_star    
    memory params.memory_star

    input:
    tuple val(sample_id), path(reads)
    path(idx)

    output:
    path("${sample_id}_Aligned.out.bam"), emit: bams
    path("${sample_id}_*.out")
    path("${sample_id}_*.tab")

    script:
    """
    STAR \
        --readFilesCommand 'gunzip -c' \
        --genomeDir $idx \
        --runThreadN $task.cpus \
        --outFileNamePrefix "${sample_id}_" \
        --outSAMtype BAM Unsorted \
        --readFilesIn $reads \
        $params.additional_star
    """

}

process FEATURECOUNTS {

    errorStrategy 'finish'

    publishDir "${params.outdir}/counts/", mode: 'copy'
    
    cpus params.cpus_featurecounts    
    memory params.memory_featurecounts

    input:
    path(bams)
    path(gtf)

    output:
    path("*.txt")
    
    script:
    """
    Rscript --vanilla $baseDir/bin/featurecounts.r $gtf $task.cpus    
    """

}

//--------------------------------------------------------------------------
// Workflow
//--------------------------------------------------------------------------

ch_fastq = Channel.fromFilePairs(params.fastq, checkIfExists: true)

workflow RNA_STAR {

    FASTQC(ch_fastq)

    if(params.run_index){
        INDEX(params.genome, params.gtf) // CI use case
        STAR(ch_fastq, INDEX.out.idx)    // CI use case
    } else {
        STAR(ch_fastq, params.index) // normal use case
    } 
    
    // if command line was --prefix_counts '' groovy will turn that into true
    if(params.prefix_counts instanceof Boolean) {
        prefix_counts = ""
    } else prefix_counts = params.prefix_counts

    FEATURECOUNTS(STAR.out.bams.collect(), params.gtf)

    // Summary message
    def od = params.outdir
    workflow.onComplete {
        Date date2 = new Date()
        String datePart2 = date2.format("yyyy-dd-MM -- ")
        String timePart2 = date2.format("HH:mm:ss")
        def end_date = datePart2 + timePart2
        ANSI_RESET   = "\u001B[0m"
        ANSI_GREEN   = "\u001B[32m"
        DASHEDDOUBLE = "=".multiply(120)

        println ""
        println "$ANSI_GREEN" + "$DASHEDDOUBLE"
        println "Pipeline completed!"
        println "End: $end_date"
        println "Results are in: " + od
        println "$DASHEDDOUBLE" + "$ANSI_RESET"
        println ""
    }


}

workflow { RNA_STAR() }