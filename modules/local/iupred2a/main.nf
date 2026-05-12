process IUPRED2A {
    tag "${meta.id}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(fasta)
    val mode
    val anchor
    path iupred2a_dir

    output:
    tuple val(meta), path('*.iupred2a.*.txt'), emit: predictions

    script:
    def anchorFlag = anchor ? '-a' : ''
    def resultLabel = anchor ? "${mode}.anchor2" : mode
    def outFile = "${meta.id}.iupred2a.${resultLabel}.txt"
    """
    python3 ${iupred2a_dir}/iupred2a.py ${anchorFlag} ${fasta} ${mode} > ${outFile}
    """
}
