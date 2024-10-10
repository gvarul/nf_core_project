process SALMON_INDEX {

    publishDir "indexed", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'biocontainers/salmon:1.10.1--h7e5ed60_0' }"


    input:
    path genome_fasta
    path read_fasta

    output:
    path "salmon", emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //argument from command prompt
    def args = task.ext.args?: ''

    //get decoy ids from genome, grep - all the lines starts with >, cut - extracts first field (seq_id), saves in decoys.txt 
    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"

    // containes genome & transcriptome sequences
    def gentrome = "gentrome.fa.gz"

    //handling gzipped fies
    if (genome_fasta.endsWith('.gz')){
        get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d '' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
        gentrome = "gentromr.fa.gz"
    }

    
    """
    $get_decoy_ids
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    //salmon Indexing

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        -d decocys.txt \\
        $args \\
        -i salmon

//version tracking 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon//g")
    END_VERSIONS
    """

    //Directory generation
    stub:
    """
    mkdir salmon
    touch salmon/complete_ref_lens.bin
    touch salmon/ctable.bin
    touch salmon/ctg_offsets.bin
    touch salmon/duplicate_clusters.tsv
    touch salmon/info.json
    touch salmon/mphf.bin
    touch salmon/pos.bin
    touch salmon/pre_indexing.log
    touch salmon/rank.bin
    touch salmon/refAccumLengths.bin
    touch salmon/ref_indexing.log
    touch salmon/reflengths.bin
    touch salmon/refseq.bin
    touch salmon/seq.bin
    touch salmon/versionInfo.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon//g")
    END_VERSIONS
    """


}

process SALMON_QUANT {
    publishDir " quantifying", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'biocontainers/salmon:1.10.1--h7e5ed60_0' }"


    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    path transcript_fasta
    val alignment_mode
    val lib_type

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val(meta), path("*info.json"), emit: json_info, optional:true
    tuple val(meta), path("*lib_format_counts.json"), emit: lib_format_counts, optional: true
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"  //output directory using meta.id

    def reference   = "--index $index"
    def reads1 = [], reads2 = []
    meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }
    def input_reads = meta.single_end ? "-r ${reads1.join(" ")}" : "-1 ${reads1.join(" ")} -2 ${reads2.join(" ")}"
    if (alignment_mode) {
        reference   = "-t $transcript_fasta"
        input_reads = "-a $reads"
    }

    def strandedness_opts = [
        'A', 'U', 'SF', 'SR',
        'IS', 'IU' , 'ISF', 'ISR',
        'OS', 'OU' , 'OSF', 'OSR',
        'MS', 'MU' , 'MSF', 'MSR'
    ]

    def strandedness =  'A'
    if (lib_type) {
        if (strandedness_opts.contains(lib_type)) {
            strandedness = lib_type
        } else {
            log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    } else {
        strandedness = meta.single_end ? 'U' : 'IU'
        if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? 'SF' : 'ISF'
        } else if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? 'SR' : 'ISR'
        }
    }

    """

    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        $reference \\
        $input_reads \\
        $args \\
        -o $prefix

    if [ -f $prefix/aux_info/meta_info.json ]; then
        cp $prefix/aux_info/meta_info.json "${prefix}_meta_info.json"
    fi
    if [ -f $prefix/lib_format_counts.json ]; then
        cp $prefix/lib_format_counts.json "${prefix}_lib_format_counts.json"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}_meta_info.json
    touch ${prefix}_lib_format_counts.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

}

workflow {
    //process flows
    SALMON_QUANT(
        meta: [id: 'sample_1'],
        reads: params.reads,
        index: params.salmon_index,
        gtf: params.gtf
        transcript_fasta: params.transcript_fasta,
        alignment_mode: params.alignment_mode,
        lib_type: params.lib_type
    )
}