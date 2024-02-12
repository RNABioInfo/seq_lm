import java.nio.file.NoSuchFileException

/**
 * Take a map of input arguments, find valid inputs, and return a channel
 * with elements of `[metamap, seqs.fastq.gz, path-to-bamstats-stats]`.
 * The last item is `null` if `bamStats` was not run.
 *
 * @param arguments: map with arguments containing
 *  - "input": path to either: (i) input BAM file, (ii) top-level directory containing
 *     BAM files, (iii) directory containing sub-directories which contain BAM
 *     files
 *  - "sample": string to name single sample
 *  - "bamstats_stats": boolean whether to write the `bamstats` stats
 * @return Channel of `[Map(alias, ...), Path, Path|null]`.
 *  The first element is a map with metadata, the second is the path to the
 *  `.bam` file with the alignments and the third is
 *  the path to the directory with the bam statistics (or `null` if `bamStats`
 *  wasn't run).
 */
def bamIngress(Map arguments) {
    // check arguments
    Map margs = parse_arguments(arguments)
    // define the channel for holding the inputs [metamap, input_path]. It will be
    // either filled by `watchPath` (only emitting files) or by the data of the three
    // input types (single file or dir with bam or subdirs with bam).
    def ch_input
    // handle `watchPath` case
    if (margs['watch_path']) {
        ch_input = watch_path(margs)
    } else {
        error 'Retrospective analysis not yet implemented.'
    // create a channel with the inputs (single file / dir with bam / subdirs
    // with bam)
    // ch_input = get_valid_inputs(margs)
    }

    def ch_result
    if (margs.bam_stats) {
        // run bamstats regardless of input type
        //ch_result = bamStats(ch_input.map { [it[0], it[1]] })
        ch_result = ch_input.map { meta, bam, allBam -> [meta, bam, allBam] }
    } else {
        // the bam stats were not requested
        ch_result = ch_input.map { meta, bam, allBam -> [meta, bam, allBam] }
    }

    return ch_result
}

/**
 * Run `watchPath` on the input directory and return a channel [metamap, path-to-bam].
 * The contains the `alias` (either `margs["sample"]` or the name of the parent
 * directory of the file).
 *
 * @param margs: map with parsed input arguments
 * @return: Channel of [metamap, path-to-bam]
 */
def watch_path(Map margs) {
    log.info "Watching path $margs.input"
    // we have one case to consider: (i) files being generated in sub-directories.
    Path input
    try {
        input = file(margs.input, checkIfExists: true)
    } catch (NoSuchFileException e) {
        error "Input path $margs.input does not exist."
    }

    if (input.isFile()) {
        error "Input ($input) must be a directory when using `watch_path`."
    }
    // get existing BAM files first (look for relevant files in the top-level dir and
    // all sub-dirs)
    def ch_existing_input = Channel.fromPath(input)
    | concat(Channel.fromPath("$input/*", type: 'dir'))
    | map { get_bam_files_in_dir(it) }
    | flatten
    | filter { it.name.endsWith('.bam') }

    // now get channel with files found by `watchPath`
    def ch_watched = Channel.watchPath("$input/**").until { it.name.startsWith('STOP') }
    // only keep BAM files
    | filter { it.name.endsWith('.bam') }
    // merge the channels
    // ch_watched = ch_existing_input | concat(ch_watched)
    // check if input is as expected; start by throwing an error when finding files in
    // top-level dir and sub-directories
    // String prev_input_type
    // ch_watched
    // | map {
    //     String input_type = (it.parent == input) ? 'top-level' : 'sub-dir'
    //     if (prev_input_type && (input_type != prev_input_type)) {
    //         error '`watchPath` found BAM files in the top-level directory ' +
    //             'as well as in sub-directories.'
    //     }
    //     // if file is in a sub-dir, make sure it's not a sub-sub-dir
    //     if ((input_type == 'sub-dir') && (it.parent.parent != input)) {
    //         error '`watchPath` found a BAM file more than one level of ' +
    //             "sub-directories deep ('$it')."
    //     }
    //     prev_input_type = input_type
    // }

    ch_watched = ch_watched
    | map {
        // This file could be in the top-level dir or a sub-dir. In the first case
        // check if a sample name was provided. In the second case, the alias is
        // always the name of the sub-dir.
        [create_metamap([runName: margs['runName'], replicateName: it.parent.name]), it, file(it.parent / '*.bam')]
    }

    return ch_watched
}

process bamStats {
    label 'preproc'
    cpus 3
    input:
        tuple val(meta), path(input)
    output:
        tuple val(meta), path(input), path('bam_stats')
    script:
        String bam_stats_outdir = 'bam_stats'
        String out = "$bam_stats_outdir/bam_stats.tsv.gz"

        """
        bamstats --help
        samtools index -@ $task.cpus $input
        mkdir $bam_stats_outdir
        bamstats \
            --histograms=$bam_stats_outdir/histograms \
            -s ${meta['alias']} \
            -f $bam_stats_outdir/flag_stats.tsv \
            -t $task.cpus \
            $input \
            | bgzip -@ $task.cpus > $out
        """
}

/**
 * Parse input arguments for `bam_ingress`.
 *
 * @param arguments: map with input arguments (see `bam_ingress` for details)
 * @return: map of parsed arguments
 */
Map parse_arguments(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args:['input'],
        kwargs:['runName': null,
                'bam_stats': false,
                'watch_path': true],
        name: 'bam_ingress')
    return parser.parse_args(arguments)
}

/**
 * Find valid inputs based on the input type.
 *
 * @param margs: parsed arguments (see `bam_ingress` for details)
 * @return: channel of `[metamap, input-path]`; `input-path` can be the path to
 *  a single BAM file or to a directory containing BAM files
 */
// def get_valid_inputs(Map margs) {
//     log.info 'Checking bam input.'
//     Path input
//     try {
//         input = file(margs.input, checkIfExists: true)
//     } catch (NoSuchFileException e) {
//         error "Input path $margs.input does not exist."
//     }
//     // declare resulting input channel and other variables needed in the outer scope
//     def ch_input
//     ArrayList sub_dirs_with_bam_files
//     // handle case of `input` being a single file
//     if (input.isFile()) {
//         // the `bamStats` process can deal with directories or single file inputs
//         ch_input = Channel.of(
//             [create_metamap([runName: margs["runName"], replicateName: replicateName]), input])
//     } else if (input.isDirectory()) {
//         // input is a directory --> we accept two cases: (i) a top-level directory with
//         // bam files and no sub-directories or (ii) a directory with one layer of
//         // sub-directories containing bam files
//         ArrayList dir_has_bam_files = get_bam_files_in_dir(input)
//         // find potential sub-directories (and sub-dirs with BAM files; note that
//         // these lists can be empty)
//         ArrayList sub_dirs = file(input.resolve('*'), type: 'dir')
//         sub_dirs_with_bam_files = sub_dirs.findAll { get_bam_files_in_dir(it) }
//         // deal with first case (top-lvl dir with BAM files and no sub-directories
//         // containing BAM files)
//         if (dir_has_bam_files) {
//             if (sub_dirs_with_bam_files) {
//                 error "Input directory '$input' cannot contain BAM " +
//                     'files and sub-directories with BAM files.'
//             }
//             ch_input = Channel.of(
//                 [create_metamap([runName: margs["runName"], replicateName: replicateName]), input])
//         } else {
//             // deal with the second case (sub-directories with bam data) --> first
//             // check whether we actually found sub-directories
//             if (!sub_dirs_with_bam_files) {
//                 error "Input directory '$input' must contain either BAM files " +
//                     'or sub-directories containing BAM files.'
//             }
//             // make sure that there are no sub-sub-directories with bam files and that
//             // the sub-directories actually contain bam files)
//             if (sub_dirs.any {
//                 ArrayList subsubdirs = file(it.resolve('*'), type: 'dir')
//                 subsubdirs.any { get_bam_files_in_dir(it) }
//             }) {
//                 error "Input directory '$input' cannot contain more " +
//                     'than one level of sub-directories with BAM files.'
//             }
//             ch_input = Channel.fromPath(sub_dirs_with_bam_files).map {
//                     [create_metamap([runName: margs["runName", replicateName: replicateName]]), it]
//             }
//         }
//     } else {
//         error "Input $input appears to be neither a file nor a directory."
//     }
//     return ch_input
// }

/**
 * Create a map that contains at least these keys: `[alias, barcode, type]`.
 * `alias` is required, `barcode` and `type` are filled with default values if
 * missing. Additional entries are allowed.
 *
 * @param kwargs: map with input parameters; must contain `alias`
 * @return: map(alias, barcode, type, ...)
 */
Map create_metamap(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args: ['runName'],
        kwargs: [
            'replicateName': null,
        ],
        name: 'create_metamap',
    )
    return parser.parse_known_args(arguments)
}

/**
 * Get the bam files in the directory (non-recursive).
 *
 * @param dir: path to the target directory
 * @return: list of found bam files
 */
ArrayList get_bam_files_in_dir(Path dir) {
    return file(dir / '*.bam', type: 'file')
}
