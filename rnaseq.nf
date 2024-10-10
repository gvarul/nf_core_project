
process TRIM_GALORE{

    publishDir "TRIMMED", mode: 'copy'

    input:
        tuple val (sampleid), path(reads)

    output:
        path "*"
        path "*trimmed*.fq.gz", emit:trimmed //to lable the files

    script:
    //Trimm anything below 20 Phred Score
    """
    trim_galore  --paired --q 20 --gzip --basename ${sampleid}_trimmed ${reads}
    """
}

process QC{
    publishDir "QC_REPORT", mode: 'copy'

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

    publishDir "REF_INDEX", mode:'copy'

    input:
        path(ref_fasta)
        path(ref_gtf)

    output:
        path "*", emit:index

    script:
    """
    STAR --runThreadN 64 \\
    --runMode genomeGenerate \\
    --genomeDir index \\
    --genomeFastaFiles ${ref_fasta} \\
    --sjdbGTFfile ${ref_gtf} \\
    --genomeSAindexNbases 12

    """
}

process STAR_MAPPING {
    publishDir "MAPPING", mode:'copy'
    cpus 96

    input:
        tuple val(sampleid), path(read1), path(read2),path(index)
    
    output:
        path "*"
        path "*.bam", emit:bams

    script:
    """
    STAR --runThreadN 72 --genomeDir ${index} \\
    --readFilesIn ${read1} ${read2} \\
    --outSAMtype BAM SortedByCoordinate \\
    --outFileNamePrefix ${sampleid} \\
    --readFilesCommand zcat
    """
}

process FEATURE_COUNT {
    publishDir "READ_COUNT", mode: 'copy'

    input:
        path(bams)
        path(ref_gtf)
        val(strand)

    output:
        path "*"

    script:
    """
    featureCounts -T 64 -s ${strand} -p --countReadPairs -t exon \\
    -g gene_id -Q 10 -a ${ref_gtf} -o gene_count ${bams}
    
    multiqc gene_count*
    """
}



workflow{

//Reference Genomes

ref_fasta=Channel.fromPath(params.ref_fasta)
ref_gtf=Channel.fromPath(params.ref_gtf)

fastq_ch=Channel.fromFilePairs(params.reads)
strand=Channel.of(params.strand)

// Trimminng 

TRIM_GALORE(fastq_ch).set{trimmed}


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

//Test
//ref_fasta.view()
//ref_gtf.view()
//strand.view()
//fastq_ch.view()
}
