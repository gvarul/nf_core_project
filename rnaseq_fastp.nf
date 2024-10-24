process TRIMMED_FASTP {

publishDir "TRIMMED", mode:'copy'

input: 
	tuple val(sampleid), path(reads)

output: 
	path "*"
	path "*trimmed*.fq.gz", emit:trimmed

script: 
"""
fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleid}_trimmed_1.fq.gz -O ${sampleid}_trimmed_2.fq.gz -q 20 -h ${sampleid}_fastp.html

"""

}


workflow {

ref_fasta=Channel.fromPath(params.ref_fasta)
ref_gtf=Channel.fromPath(params.ref_gtf)

fastq_ch=Channel.fromFilePairs(params.reads)

strand=Channel.of(params.strand)



// Trimminng 
TRIMMED_FASTP(fastq_ch).set{trimmed}

}



//Completed at: 10-Oct-2024 14:14:34
//Duration    : 1h 16m 41s
//CPU hours   : 8.3
//Succeeded   : 12