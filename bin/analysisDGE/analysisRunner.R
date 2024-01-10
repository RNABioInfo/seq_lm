library(optparse)

# Create an option parser
parser <- OptionParser()
# Add options for bamPath, annotationFile, and outputPath
parser <- add_option(parser, "-b", "--bamPath", dest = "bamPath", help = "Path to the BAM file")
parser <- add_option(parser, "-a", "--annotationFile", dest = "annotationFile", help = "Path to the annotation GTF/GFF file")
parser <- add_option(parser, "-o", "--outputPath", dest = "outputPath", help = "Path to the output file")

# Parse the command line arguments
options <- parse_args(parser)

# Access the values of the options
bamPath <- options$bamPath
annotationFile <- options$annotationFile
outputPath <- options$outputPath

# Use the values in your code
# ...
