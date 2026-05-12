nextflow.enable.dsl = 2

include { IUPRED2A } from './modules/local/iupred2a/main'

def parseBooleanParam(value) {
    if (value instanceof Boolean) {
        return value
    }

    def normalized = value?.toString()?.trim()?.toLowerCase()
    if (['true', '1', 'yes', 'y'].contains(normalized)) {
        return true
    }
    if (['false', '0', 'no', 'n', ''].contains(normalized)) {
        return false
    }

    throw new IllegalArgumentException("Invalid boolean value for --anchor: ${value}")
}

workflow {
    if (!params.input) {
        error "Missing required parameter --input. Example: --input 'iupred2a/P53_HUMAN.seq'"
    }

    def validModes = ['long', 'short', 'glob']
    def mode = params.mode.toString().trim().toLowerCase()
    if (!validModes.contains(mode)) {
        error "Invalid --mode '${params.mode}'. Expected one of: ${validModes.join(', ')}"
    }

    def anchor = parseBooleanParam(params.anchor)
    def iupred2aDir = file("${projectDir}/iupred2a")

    ch_fasta = channel
        .fromPath(params.input, checkIfExists: true)
        .map { fasta ->
            def sampleId = fasta.baseName.replaceAll(/[^A-Za-z0-9_.-]/, '_')
            tuple([id: sampleId], fasta)
        }

    IUPRED2A(ch_fasta, mode, anchor, iupred2aDir)

    IUPRED2A.out.predictions.view { meta, prediction ->
        "IUPred2A result for ${meta.id}: ${workflow.launchDir}/${params.outdir}/${prediction.name}"
    }
}
