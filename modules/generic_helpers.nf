nextflow.enable.types = true

process getVersions {
    label 'preproc'
    cpus 1
    output:
        file('versions.txt')
    script:
        """
        python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
        """
}

process getParams {
    label 'seqLM'
    cpus 1
    output:
        file('params.json')
    script:
        def paramsJSON = new groovy.json.JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    debug true
    label 'seqLM'
    publishDir(
        params.ex_dir,
        mode: 'copy',
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple(fname: Path, dirname: String?)
    output:
        file(fname.name)
    script:
        '''
        '''
}
