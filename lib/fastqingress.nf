process handleSingleFile {
    label "microbial"
    cpus 1
    input:
        file reads
    output:
        path "$reads.simpleName"
    script:
        def name = reads.simpleName
        def reads_dir = 'reads_dir'
    """
    mkdir $name
    mv $reads $name
    """
}


process checkSampleSheet {
    label "microbial"
    cpus 1
    input:
        file "sample_sheet.txt"
    output:
        file "samples.txt"
    """
    check_sample_sheet.py sample_sheet.txt samples.txt
    """
}


/**
 * Take an input file and sample name to return a channel with
 * a single named sample.
 *
 *
 * @param input_file Single fastq file
 * @param sample_name Name to give the sample
 * @return Channel of tuples (path, sample_id, type)
 */
def handle_single_file(input_file, sample_name)
{
    singleFile = Channel.fromPath(input_file)
    sample = handleSingleFile(singleFile)
    return sample.map { it -> tuple(it, sample_name === null ? it.simpleName : sample_name, 'test_sample') }
}


/**
 * Find fastq data using various globs. Wrapper around Nextflow `file`
 * method.
 *
 * @param pattern file object corresponding to top level input folder.
 * @param maxdepth maximum depth to traverse
 * @return list of files.
 */

def find_fastq(pattern, maxdepth)
{
    files = []
    extensions = ["fastq", "fastq.gz", "fq", "fq.gz"]
    for (ext in extensions) {
        files += file(pattern.resolve("*.${ext}"), type: 'file', maxdepth: maxdepth)
    }
    return files
}


/**
 * Rework EPI2ME flattened directory structure into standard form
 * files are matched on barcode\d+ and moved into corresponding
 * subdirectories ready for processing.
 *
 * @param input_folder Top-level input directory.
 * @param staging Top-level output_directory.
 * @return A File object representating the staging directory created
 *     under output_folder
 */
def sanitize_fastq(input_folder, staging)
{
    // TODO: this fails if input_folder is an S3 path
    println("Running sanitization.")
    println(" - Moving files: ${input_folder} -> ${staging}")
    staging.mkdirs()
    files = find_fastq(input_folder.resolve("**"), 1)
    for (fastq in files) {
        fname = fastq.getFileName()
        // find barcode
        pattern = ~/barcode\d+/
        matcher = fname =~ pattern
        if (!matcher.find()) {
            // not barcoded - leave alone
            fastq.renameTo(staging.resolve(fname))
        } else {
            bc_dir = file(staging.resolve(matcher[0]))
            bc_dir.mkdirs()
            fastq.renameTo(staging.resolve("${matcher[0]}/${fname}"))
        }
    }
    println(" - Finished sanitization.")
    return staging
}


/**
 * Take an input directory return the barcode and non barcode
 * sub directories contained within.
 *
 *
 * @param input_directory Top level input folder to locate sub directories
 * @return A list containing sublists of barcode and non_barcode sub directories
 */
def get_subdirectories(input_directory)
{
    barcode_dirs = file("$input_directory/barcode*", type: 'dir', maxdepth: 1)
    all_dirs = file("$input_directory/*", type: 'dir', maxdepth: 1)
    non_barcoded = ( all_dirs + barcode_dirs ) - all_dirs.intersect(barcode_dirs)
    return [barcode_dirs, non_barcoded]
}


/**
 * Load a sample sheet into a Nextflow channel to map barcodes
 * to sample names.
 *
 * @param samples CSV file according to MinKNOW sample sheet specification
 * @return A Nextflow Channel of tuples (barcode, sample name, sample type)
 */
def get_sample_sheet(sample_sheet)
{
    println("Checking sample sheet.")
    sample_sheet = file(sample_sheet);
    is_file = sample_sheet.isFile()

    if (!is_file) {
        println('Error: `--samples` is not a file.')
        exit 1
    }

    return checkSampleSheet(sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.barcode,
            row.sample_id,
            row.type ? row.type : 'test_sample')
        }
}


/**
 * Take a list of input directories and return directories which are
 * valid, i.e. contains only .fastq(.gz) files.
 *
 *
 * @param input_dirs List of barcoded directories (barcodeXX,mydir...)
 * @return List of valid directories
 */
def get_valid_directories(input_dirs)
{
    valid_dirs = []
    no_fastq_dirs = []
    invalid_files_dirs = []
    for (d in input_dirs) {
        valid = true
        fastq = find_fastq(d, 1)
        all_files = file(d.resolve("*"), type: 'file', maxdepth: 1)
        non_fastq = ( all_files + fastq ) - all_files.intersect(fastq)

        if (non_fastq) {
            valid = false
            invalid_files_dirs << d
        }
        if (!fastq) {
            valid = false
            no_fastq_dirs << d
        }
        if (valid) {
            valid_dirs << d
        }
    }
    if (valid_dirs.size() == 0) {
        error_message = "Error: None of the directories given contain .fastq(.gz) files."
        println(error_message)
        exit 1
    }
    if (no_fastq_dirs.size() > 0) {
        println("Warning: Excluding directories not containing .fastq(.gz) files:")
        for (d in no_fastq_dirs) {
            println("   - ${d}")
        }
    }
    if (invalid_files_dirs.size() > 0) {
        println("Warning: Excluding directories containing non .fastq(.gz) files:")
        for (d in invalid_files_dirs) {
            println("   - ${d}")
        }
    }
    return valid_dirs
}


/**
 * Take an input directory and sample name to return a channel
 * with a single named sample.
 *
 *
 * @param input_directory Directory of fastq files
 * @param sample_name Name to give the sample
 * @return Channel of tuples (path, sample_id, type)
 */
def handle_flat_dir(input_directory, sample_name)
{
    valid_dirs= get_valid_directories([ file(input_directory) ])
    return Channel.fromPath(valid_dirs)
        .map { it -> tuple(it, sample_name === null ? it.baseName : sample_name, 'test_sample') }
}


/**
 * Take a list of barcode directories and a sample sheet to return
 * a channel of named samples.
 *
 *
 * @param barcoded_dirs List of barcoded directories (barcodeXX,...)
 * @param sample_sheet List of tuples mapping barcode to sample name
 *     or a simple string for non-multiplexed data.
 * @return Channel of tuples (path, sample_id, type)
 */
def handle_barcoded_dirs(barcoded_dirs, sample_sheet)
{
    valid_dirs = get_valid_directories(barcoded_dirs)
    // link sample names to barcode through sample sheet
    if (!sample_sheet) {
        sample_sheet = Channel
            .fromPath(valid_dirs)
            .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
            .map { path -> tuple('test_sample', '', path.baseName, path.baseName) }
    }
    return Channel
        .fromPath(valid_dirs)
        .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
        .map { path -> tuple(path.baseName, path) }
        .join(sample_sheet)
        .map { barcode, path, sample, type -> tuple(sample, barcode, path, type) }
}


/**
 * Take a list of non-barcode directories to return a channel
 * of named samples. Samples are named by directory baseName.
 *
 *
 * @param non_barcoded_dirs List of directories (mydir,...)
 * @return Channel of tuples (path, sample_id, type)
 */
def handle_non_barcoded_dirs(non_barcoded_dirs)
{
    valid_dirs = get_valid_directories(non_barcoded_dirs)
    return Channel.fromPath(valid_dirs)
        .map { path -> tuple(path.baseName, path, path.baseName, 'test_sample') }
}


/**
 * Take an input (file or directory) and return a channel of
 * named samples.
 *
 *
 * @param input Top level input file or folder to locate fastq data.
 * @param sample string to name single sample data.
 * @param sample_sheet Path to sample sheet CSV file.
 * @return Channel of tuples (path, sample_id, type)
 */
def fastq_ingress(input, output_folder, sample, sample_sheet, sanitize)
{
    println("Checking fastq input.")
    input = file(input);

    // Handle file input
    if (input.isFile()) {
        // Assume sample is a string at this point
        println('Single file input detected.')
        if (sample_sheet) {
            println('Warning: `--sample_sheet` given but single file input found. Ignoring.')
        }
        return handle_single_file(input, sample)
    }

    // Handle directory input
    if (input.isDirectory()) {
        // EPI2ME harness
        if (sanitize) {
            staging = file(output_folder).resolve("staging")
            input = sanitize_fastq(file(input), staging)
        }

        // Get barcoded and non barcoded subdirectories
        (barcoded, non_barcoded) = get_subdirectories(input)

        // Case 03: If no subdirectories, handle the single dir
        if (!barcoded && !non_barcoded) {
            println("Single directory input detected.")
            if (sample_sheet) {
                println('Warning: `--sample_sheet` given but single non-barcode directory found. Ignoring.')
            }
            return handle_flat_dir(input, sample)
        }

        if (sample) {
            println('Warning: `--sample` given but multiple directories found, ignoring.')
        }

        // Case 01, 02, 04: Handle barcoded and non_barcoded dirs
        // Handle barcoded folders
        barcoded_samples = Channel.empty()
        if (barcoded) {
            println("Barcoded directories detected.")
            if (sample_sheet) {
                sample_sheet = get_sample_sheet(sample_sheet)
            }
            barcoded_samples = handle_barcoded_dirs(barcoded, sample_sheet)
        }

        non_barcoded_samples = Channel.empty()
        if (non_barcoded) {
            println("Non barcoded directories detected.")
            if (!barcoded && sample_sheet) {
                println('Warning: `--sample_sheet` given but no barcode directories found.')
            }
            non_barcoded_samples = handle_non_barcoded_dirs(non_barcoded)
        }

        return barcoded_samples.mix(non_barcoded_samples)
    }
}
