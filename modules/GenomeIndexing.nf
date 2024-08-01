process runGenomeIndexing {
    publishDir "${params.index_dir}", mode: 'move'
    
    output:
    val true  // for state dependency
    path "*"

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles $params.fasta_file \
    --sjdbGTFfile $params.annotation_file \
    --limitGenomeGenerateRAM $params.max_RAM_indexing \
    --runThreadN $params.genome_indexing_nt
    """
}
