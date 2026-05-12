nextflow.enable.dsl = 2

include { IUPRED2A } from './modules/local/iupred2a/main'

process SPLIT_FASTA {
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)

    output:
    path '*.fa', emit: records

    script:
    def sampleId = meta.id
    """
    mkdir -p split_records

    awk '
    function sanitize(raw, fallback, value) {
        value = raw
        sub(/^[[:space:]]+/, "", value)
        sub(/[[:space:]].*\$/, "", value)
        gsub(/[^A-Za-z0-9_.-]/, "_", value)
        gsub(/^_+/, "", value)
        gsub(/_+\$/, "", value)
        if (value == "") {
            value = fallback
        }
        return value
    }

    function close_record() {
        if (out_file != "") {
            close(out_file)
        }
    }

    /^>/ {
        close_record()
        record_count++
        sequence_id = sanitize(substr(\$0, 2), sprintf("seq%06d", record_count))
        out_file = sprintf("split_records/%06d__%s.fa", record_count, sequence_id)
        print ">" sequence_id > out_file
        next
    }

    NF {
        if (out_file == "") {
            record_count++
            sequence_id = sprintf("seq%06d", record_count)
            out_file = sprintf("split_records/%06d__%s.fa", record_count, sequence_id)
            print ">" sequence_id > out_file
        }
        gsub(/[[:space:]]/, "", \$0)
        if (\$0 != "") {
            print toupper(\$0) >> out_file
        }
    }

    END {
        close_record()
        if (record_count == 0) {
            exit 1
        }
    }
    ' ${fasta}

    record_count=\$(find split_records -maxdepth 1 -type f -name '*.fa' | wc -l)
    if [ "\${record_count}" -eq 1 ]; then
        cp split_records/*.fa "${sampleId}.fa"
    else
        for record in split_records/*.fa; do
            cp "\${record}" "${sampleId}__\$(basename "\${record}")"
        done
    fi
    """
}

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
    def publishedOutdir = params.outdir.toString()
    if (!publishedOutdir.startsWith('/')) {
        publishedOutdir = "${workflow.launchDir}/${publishedOutdir}"
    }

    ch_fasta_files = channel
        .fromPath(params.input, checkIfExists: true)
        .map { fasta ->
            def sampleId = fasta.baseName.replaceAll(/[^A-Za-z0-9_.-]/, '_')
            tuple([id: sampleId], fasta)
        }

    SPLIT_FASTA(ch_fasta_files)

    ch_fasta_records = SPLIT_FASTA.out.records
        .flatten()
        .map { fasta ->
            tuple([id: fasta.baseName], fasta)
        }

    IUPRED2A(ch_fasta_records, mode, anchor, iupred2aDir)

    IUPRED2A.out.predictions.view { meta, prediction ->
        "IUPred2A result for ${meta.id}: ${publishedOutdir}/${prediction.name}"
    }
}
