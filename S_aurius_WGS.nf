#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/****************************************************
 * Parameters
 ****************************************************/
params.accessions = [
    "SRR27321999","SRR27322000","SRR27322001","SRR27322002",
    "SRR27322003","SRR27322004","SRR27322005","SRR27322006",
    "SRR27322007","SRR27322008"
]

params.outdir = "./results"

/****************************************************
 * Processes
 ****************************************************/

/* Step 1: Download + convert SRA → FASTQ */
process sra_download {
    tag "$sra_id"
    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}*.fastq.gz")

    script:
    """
    set -euo pipefail
    module load sratoolkit/3.0.0
    echo "Downloading: ${sra_id}"
    prefetch ${sra_id}
    fasterq-dump --outdir . --split-3 ${sra_id}/${sra_id}.sra
    gzip ${sra_id}*.fastq
    """
}

/* Step 2: FastQC */
process fastqc {
    tag { sample_id }
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip

    script:
    """
    module load fastqc
    fastqc -t ${task.cpus} -o . ${reads}
    """
}

process multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path fastqc_reports

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc ./ -n multiqc_report.html
    """
}

/* Step 3: Trimming */
process trim_galore {
    tag { sample_id }
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz")

    script:
    """
    trim_galore --paired ${reads[0]} ${reads[1]} -o .
    """
}

/* Step 4: Assembly */
process shovill {
    tag { sample_id }
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("shovill_output/contigs.fa")

    script:
    """
    shovill --R1 ${reads[0]} --R2 ${reads[1]} --outdir shovill_output --cpus ${task.cpus}
    """
}

/* Step 5: QC of assemblies */
process quast {
    tag { sample_id }
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("quast_results/report.html")

    script:
    """
    quast.py -o quast_results ${assembly}
    """
}

process busco {
    tag { sample_id }
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("busco_output/short_summary.*.txt")

    script:
    """
    busco -i ${assembly} -o busco_output -m genome -l ${params.busco_db} -c ${task.cpus}
    """
}

/* Step 6: Annotation + Typing */
process prokka {
    tag { sample_id }
    publishDir "${params.outdir}/prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("prokka_output")

	// activate conda environment
    conda '/home/user/anaconda3/envs/prokka'
	
    script:
    """
    prokka --outdir prokka_output --prefix ${sample_id} ${assembly}
    """
}

process mlst {
    tag { sample_id }
    publishDir "${params.outdir}/mlst", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("*.mlst")

	// activate conda environment
    conda '/home/user/anaconda3/envs/mlst'
	
    script:
    """
    mlst ${assembly} > ${sample_id}.mlst
    """
}

process spatyper {
    tag { sample_id }
    publishDir "${params.outdir}/spatyper", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("*.spa")

    script:
    """
    spatyper -i ${assembly} -o ${sample_id}.spa
    """
}

process agrvate {
    tag { sample_id }
    publishDir "${params.outdir}/agrvate", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("*.agr")

	// activate conda environment
    conda '/home/user/anaconda3/envs/agrvate'
	
    script:
    """
    agrvate -i ${assembly} -o ${sample_id}.agr
    """
}

/* Step 7: AMR + Virulence */
process abricate {
    tag { sample_id }
    publishDir "${params.outdir}/abricate", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("*.tab")

	// activate conda environment
	conda '/home/user/anaconda3/envs/abricate'
	
    script:
    """
    abricate --db ncbi ${assembly} > ${sample_id}_resistance.tab
    abricate --db plasmidfinder ${assembly} > ${sample_id}_plasmids.tab
    """
}

process resfinder {
    tag { sample_id }
    publishDir "${params.outdir}/resfinder", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("*.resfinder")

    script:
    """
    resfinder.py -ifa ${assembly} -o . -s enterobacterales -l 0.6 -t 0.8
    """
}

process mob_suite {
    tag { sample_id }
    publishDir "${params.outdir}/mob_suite", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("mobtyper_results")

	// activate conda environment
	conda '/home/user/anaconda3/envs/abricate'
	
    script:
    """
    mob_typer -i ${assembly} -o mobtyper_results
    """
}

process phispy {
    tag { sample_id }
    publishDir "${params.outdir}/phispy", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("phispy_results")

    script:
    """
    phispy ${assembly} phispy_results
    """
}

/****************************************************
 * Workflow
 ****************************************************/

// Workflow definition
workflow {

    // Input: create channel of 10 SRA accessions
    samples = Channel.fromList(params.accessions)

    // Step 1: Download SRA → FASTQ
    reads = samples | sra_download

    // Step 2: QC
    fastqc_out = fastqc(reads)
    all_reports = fastqc_out.html.mix(fastqc_out.zip).collect()
    multiqc(all_reports)

    // Step 3: Trimming
    trimmed_reads = trim_galore(reads)

    // Step 4: Assembly
    assemblies = shovill(trimmed_reads)

    // Step 5: Assembly QC
    quast(assemblies)
    busco(assemblies)

    // Step 6: Annotation + typing
    prokka(assemblies)
    mlst(assemblies)
    spatyper(assemblies)
    agrvate(assemblies)

    // Step 7: AMR + virulence + plasmids + phages
    abricate(assemblies)
    resfinder(assemblies)
    mob_suite(assemblies)
    phispy(assemblies)
}


