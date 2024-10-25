process TRIM_GALORE {

publishDir "TRIMMED", mode:'copy'

input: 
	tuple val(sampleid), path(reads)

output: 
	path "*"
	path "*trimmed*.fq.gz", emit:trimmed

script: 
"""
trim_galore --paired -q 20 --gzip --basename ${sampleid}_trimmed ${reads}

"""

}







process QC {

publishDir "QC_REPORT", mode:'copy'

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

process SALMON_INDEX {

publishDir "INDEX", mode:'copy'

input: 
	path(fasta)// Reference transcriptome fa

output: 
	path "*", emit:index

script:
"""
salmon index -t ${fasta} -i salmon_index -k 31

"""
//--gencode Gencode transcriptome (human or mouse data with Gencode annotations)


}


process SALMON_QUANT {

publishDir "QUANT", mode:'copy'


input: 
	tuple val(sampleid),path(read1), path(read2), path(index)
	
output: 
	path "*", emit:bams
	
script:
"""
salmon quant -i ${index} -l A \\
    -1 ${read1} -2 ${read2} \\
    -p 4 --gcBias --seqBias --validateMappings \\
    -o ${sampleid}_quant
"""



}







workflow {

ref_fasta=Channel.fromPath(params.ref_fasta)
ref_gtf=Channel.fromPath(params.ref_gtf)

fastq_ch=Channel.fromFilePairs(params.reads)

strand=Channel.of(params.strand)



// Trimminng 
TRIM_GALORE(fastq_ch).set{trimmed}


// FASRQC & MULTIQC
raw_fastq=fastq_ch.map{items -> items[1] }.flatten().collect()

trimmed_fastq=trimmed.trimmed.flatten().collect()

raw_fastq.mix(trimmed_fastq).collect() | QC



// Build Salmon index
SALMON_INDEX(ref_fasta).set {salmon_index }


// Salmon quantification using trimmed reads
trimmed.trimmed.map {read1, read2 -> tuple("${read1.getFileName()}".split("_trimmed")[0], read1, read2) }
| combine(salmon_index.index) | SALMON_QUANT

}

