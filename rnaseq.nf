
def samplesheet = file(params.sample_sheet)

    Channel
        .from(samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sampleid, file(row.read1), file(row.read2)) }
        .set { fastq_ch }

process VERSIONS {
    publishDir "OUTPUT", mode: 'copy'

    container 'container community.wave.seqera.io/library/fastp_fastqc_multiqc_star_pruned:040af6e70566a925'

    output:
    path "versions.txt"

    script:
    """
    echo "Tool & Versions " > versions.txt
    echo "=========================" >> versions.txt
    echo "" >> versions.txt

    echo "multiqc Version:" >> versions.txt
    fastqc --version >> versions.txt
    echo "" >> versions.txt

    echo "fastqc Version:" >> versions.txt
    multiqc --version >> versions.txt
    echo "" >> versions.txt

    echo "trim_galore Version:" >> versions.txt
    trim_galore --version >> versions.txt
    echo "" >> versions.txt

    echo "fastp Version:" >> versions.txt
    fastp --version >> versions.txt
    echo "" >> versions.txt

    echo "STAR Version:" >> versions.txt
    STAR --version >> versions.txt
    echo "" >> versions.txt

    echo "featureCounts Version:" >> versions.txt
    featureCounts -v >> versions.txt
    echo "" >> versions.txt
    


    """
}

process TRIM_GALORE {
    publishDir "OUTPUT/TRIMMED", mode: 'copy'

    container 'container community.wave.seqera.io/library/fastp_fastqc_multiqc_star_pruned:040af6e70566a925'


    input:
        tuple val(sampleid), path(read1), path(read2)

    output:
        path "*trimmed*.fq.gz", emit: trimmed // Emit the trimmed reads

    script:
    // Trimm anything below 20 Phred Score
    """
    trim_galore --paired --q 20 --gzip --basename ${sampleid}_trimmed ${read1} ${read2}
    """
}

process TRIM_FASTP {

container 'container community.wave.seqera.io/library/fastp_fastqc_multiqc_star_pruned:040af6e70566a925'

publishDir "OUTPUT/TRIMMED", mode:'copy'

input: 
	tuple val(sampleid), path(read1), path(read2)

output: 
	path "*"
	path "*trimmed*.fq.gz", emit:trimmed

script: 
"""
fastp -i ${read1} -I ${read2} -o ${sampleid}_trimmed_1.fq.gz -O ${sampleid}_trimmed_2.fq.gz -q 20 -h ${sampleid}_fastp.html

"""

}

process QC{
    publishDir "OUTPUT/QC_REPORT", mode: 'copy'

    container 'container community.wave.seqera.io/library/fastp_fastqc_multiqc_star_pruned:040af6e70566a925'


    input:
        path(reads)

    output:
        path "*"

    script:
    """
    fastqc ${reads}
    multiqc *fastqc*

    mkdir FASTQC
    mv *fastqc* FASTQC

    """
}

process STAR_INDEX {

    publishDir "OUTPUT/REF_INDEX", mode:'copy'

    container 'container community.wave.seqera.io/library/fastp_fastqc_multiqc_star_pruned:040af6e70566a925'

    input:
        path(ref_fasta)
        path(ref_gtf)

    output:
        path "*", emit:index

    script:
    """
    STAR --runThreadN 3 \\
    --runMode genomeGenerate \\
    --genomeDir index \\
    --genomeFastaFiles ${ref_fasta} \\
    --sjdbGTFfile ${ref_gtf} \\
    --genomeSAindexNbases 12
    """
}

process STAR_MAPPING {
    publishDir "OUTPUT/MAPPING", mode: 'copy'
    
    container 'container community.wave.seqera.io/library/fastp_fastqc_multiqc_star_pruned:040af6e70566a925'

    cpus 3

    input:
        tuple val(sampleid), path(read1), path(read2), path(index)

    output:
        path "*.bam", emit: bams

    script:
    """
    STAR --runThreadN 3 --genomeDir ${index} \\
    --readFilesIn ${read1} ${read2} \\
    --outSAMtype BAM SortedByCoordinate \\
    --outFileNamePrefix ${sampleid} \\
    --readFilesCommand zcat
    """
}

process FEATURE_COUNT {
    publishDir "OUTPUT/READ_COUNT", mode: 'copy'

    container 'container community.wave.seqera.io/library/fastp_fastqc_multiqc_star_pruned:040af6e70566a925'


    input:
        path(bams)
        path(ref_gtf)
        val(strand)

    output:
        path "*"

    script:
    """
    featureCounts -T 3 -s ${strand} -p --countReadPairs -t exon \\
    -g gene_id -Q 10 -a ${ref_gtf} -o gene_count ${bams}
    
    multiqc gene_count*
    """
}


workflow {

    // Reference Genomes
    ref_fasta = Channel.fromPath(params.ref_fasta)
    ref_gtf = Channel.fromPath(params.ref_gtf)

    //fastq_ch = Channel.fromFilePairs(params.reads)
    strand = Channel.of(params.strand)

    VERSIONS()

// Trimming
    if (params.use_fastp) {
        TRIM_FASTP(fastq_ch).set { trimmed }
    } else {
        TRIM_GALORE(fastq_ch).set { trimmed }
    }


// FASRQC & MULTIQC

    raw_fastq=fastq_ch.map{it -> it[1]}.flatten().collect()
    trimmed_fastq=trimmed.trimmed.flatten().collect()
    raw_fastq.mix(trimmed_fastq).collect() | QC 





    STAR_INDEX(ref_fasta,ref_gtf).set{star_index}

    trimmed.trimmed.map{read1,read2 -> tuple("${read1.getFileName()}".split("_trimmed")[0],read1,read2 ) }
    |combine(star_index.index)|STAR_MAPPING
    |set{bams}


    bams.bams.collect().set{finalbams}

    FEATURE_COUNT (finalbams,ref_gtf,strand)

}
