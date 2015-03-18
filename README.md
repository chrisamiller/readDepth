An R package for inferring copy number changes from sequencing data

## Description

The readDepth package for R can detect copy number aberrations by measuring the depth of coverage obtained by massively parallel sequencing of the genome. It achieves higher accuracy than many other packages, and runs faster by utilizing multi-core architectures to parallelize the processing of these large data sets.

In contrast to other published methods, readDepth does not require the sequencing of a reference sample, and uses a robust statistical model that accounts for overdispersed data. It includes a method for effectively increasing the resolution obtained from low-coverage experiments by utilizing breakpoint information from paired end sequencing to do positional refinement. It can also be used to infer copy number using reads obtained from bisulfite sequencing experiments.

For a full description of the method and applications, see:
Miller, CA, et al. [ReadDepth: A Parallel R Package for Detecting Copy Number Alterations from Short Sequencing Reads](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016327). PLoS One. doi:10.1371/journal.pone.001632

#### :exclamation: Notice 
ReadDepth continues to have niche uses, especially on model organisms, but has largely been made obsolete on human data. For unmatched human samples, I now recommend using another method, such as CNVator. For calling somatic CN events from matched tumor/normal pairs, I recommend [copyCat](http://github.com/chrisamiller/copyCat).

## Contents
- [Installation](#installation)
- [Directory Setup](#directory_setup)
- [Annotations](#annotations)
- [Parameters](#parameters)
- [Usage](#usage)
- [FAQ](#faq)

## <a name="installation"></a>Installation instructions:

    #install a few packages from bioconductor
    source("http://bioconductor.org/biocLite.R")
    biocLite(c("IRanges","foreach","doMC","DNAcopy"))
    #install devtools if you don't have it already
    install.packages("devtools")
    library(devtools)
    install_github("chrisamiller/readDepth")

If you prefer to build the package by hand, follow these steps:

- Make sure that you have the dependencies from the CRAN and BioConductor repos:

- Download and build from source:

        git clone git@github.com:chrisamiller/readDepth.git
        cd readDepth/
        R CMD build readDepth
        R CMD INSTALL readDepth_0.9.8.4.tar.gz


### <a name="directory_setup"></a>Directory Setup
Start by creating a directory to hold all of your data and results. Within it, readDepth requires three sub-directories:
- reads/
This will contain all of your mapped reads, in one-based bed format, broken into one file for each chromosome.
- output/
This will be initially empty, and the tool will place the results of your analysis here
- annotations/
This will contain annotation files required by readDepth, including a list of entrypoints, and information on GC-content and mapability for the read length that you're working with. See the "Annotations" section below.
- The resulting directory structure will look like this:

    <working directory>
    +- annotations/
    |  +- gcWinds/
    |  |  +- chr1.gc.gz
    |  |  +- chr2.gc.gz
    |  |  +- . . .
    |  |
    |  +- mapability/
    |  |  +- chr1.dat.gz
    |  |  +- chr2.dat.gz
    |  |  +- . . .
    |  |
    |  +- entrypoints
    |
    +- output/
    |
    +- params
    |
    +- reads/
    |  +- chr1.bed
    |  +- chr2.bed
    |  . . .

### <a name="annotations"></a>Annotations
Annotations for common read lengths have been pre-computed for reference genomes hg18 and hg19. They can be downloaded and placed into the appropriate spot using the getAnnotations() function. Alternately, they can be manually accessed from the downloads page, copied to the annotations/ directory and untarred.

Instructions on computing additional annotations for additional read lengths can be found on the Annotations page

### <a name="parameters"></a>Parameters
The main analysis directory should also contain a tab-delimited file named "params". In it will be some or all of the following options:

    readLength      - the mean read length. (required)
   
    fdr             - the fdr rate to use.  A value of 0.01 should
                      be a good tradeoff between sensitivity and specificity
                      for most applications. (required)
 

    overDispersion  - the amount of overdispersion seen in the read
                      distribution. For Illumina reads from GAI and II machines,
                      a value of 3 works well (required). A value of 1 is 
                      equivalent to no overdispersion, or a Poisson distribution.
 

    gcWindowSize    - The smallest window size used for GC correction and the 
                      minimum window size and. If you change this value, you'll 
                      need to create new annotation files, since the provided 
                      files use a window size of 100.
 

    percCNGain      - the estimated amount of copy-number gain in this genome.
                      If you have some knowledge of the genome and, for example,
                      know that 25% of the genome is duplicated, then setting
                      this value to 0.25 will produce a more accurate model
                      and better calls. If unknown, then a value of 0.05
                      produces good results. (required)
 

    percCNLoss      - the estimated amount of copy-number gain in this genome.
                      (see percCNGain) (required)
 

    chunkSize       - the number of reads to process at one time. A larger chunk
                      size results in faster execution, but increases memory
                      usage. If you have at least 2GB of RAM per core, then we
                      suggest a value of 5e6. If you'll have 4GB per core, try
                      a value of 1e7 (required)
 

    maxCores        - the number of cores to use. If this parameter is not
                      present, the machine will use the maximum number of
                      cores available to it. (optional)
 

    readCores       - In some cases, disk IO may be the limiting resource and
                      it may be advantageous to set the number of cores used
                      while reading data lower than the number present in the
                      machine. The maxCores value will still be used for
                      operations that are not disk-intensive. (optional)
 

    verbose        -  (TRUE/FALSE) print more output to screen and plots to output
                      directory (optional - default FALSE)
 

You can also download an [example params file](http://www.example.com/).

### <a name="usage"></a>Usage
Start R, then run something like the following set of commands:

    # load the library
    library("readDepth")
    
    # create a readDepth object, then fill it by
    # reading in the params, setting up the environment,
    # creating the model, and choosing optimal bin size
    rdo = new("rdObject")
    
    # calculate depth of coverage in each bin
    rdo = readDepth(rdo)
    
    # correct the reads for mapability. This example uses a conservative
    # threshold of 0.75. In other words, if a bin is less than 75% mapable,
    # it's depth is set to NA. This prevents overcorrection.
    rdo.mapCor = rd.mapCorrect(rdo, minMapability=0.75)
    
    # do LOESS-based GC correction.
    rdo.mapCor.gcCor = rd.gcCorrect(rdo.mapCor)
    
    # segment the data using CBS. If you notice artifacts in the output, such
    # as regions of gain that span centromeres, you might try adding the
    # "rmGaps=FALSE" parameter. If you're using data with very high coverage
    # (say, greater than 10x), consider ad

## <a name="faq"></a>FAQs

- I have paired end data and would like to refine segment edges based on breakpoints called from my data. How do I generate the breakpoints.dat file referred to in the documentation?

The breakpoints.dat file is a simple 3-column bed file, listing the locations of breakpoints. readDepth does not call breakpoints - for that, you'll have to use another tool. One I use frequently is BreakDancer, but there are others as well.

- Can you provide annotation files for #bp read lengths on genome build ##?

If your avg read length is within a few base pairs of one of the existing annotation sets, you can just use that - the differences will be minor. If your data is wildly different, then check out the Annotations page for instructions and scripts that will help you create your own.

- How do I generate my own annotation files, perhaps for non-human genomes?

Check out the Annotations section above for more info on this. In short, the mapability files are created by creating all possible reads from the genome, then mapping them back using BWA. This information is then used to calculate the mapability and GC content of all mapable reads in 100 bp windows. Each read's position is specified by it's left-most base.

- Do my input bed files need to be sorted?

Yes.

- What are the columns in the output files (segs.dat or alts.dat)

The columns are chr, start, stop, number of bins, absolute copy number call

- For mate-pair or paired-end data, should I use both ends of these reads?

Since the positions of the two ends is dependent (except in the case of sv), adding the second read doesn't really give you much information. I recommend only using one end of paired reads (along with any singlet reads).

- What's the difference between the output files 'segs.dat' and 'alts.dat'?

The segs.dat file contains segmented copy number calls for all regions of the genome. The alts.dat file contains the subset of those segments that exceed the gain and loss thresholds. Also take note of the next item:

- When I run readDepth on male samples, why does it always appear that the X and Y chromosomes are entirely deleted?

In the alts.dat file, readDepth currently outputs all segments that differ from an expected ploidy of two. It's easy enough to filter this data post-hoc. We plan to handle this more intelligently in the future.

- I'm not happy with the resolution of my output segments. Is there any way to get higher resolution?

The easiest way to do that would be to change the FDR setting in the params file. If you're willing to tolerate a higher false discovery rate, you can get smaller bins. The default is usually 0.01, try 0.05.

The other thing you should probably do is check the pdfs in your output directory (assuming you set verbose to TRUE). They give an indication of how well the model fits the data. (dark bars = model, white bars = data). If you're finding that your data has less than the default 3x overdispersion, you may be able to lower the overdispersion parameter. This should also increase resolution.

- Can I use readDepth for exome data?

No. [This answer explains why](https://www.biostars.org/p/17820/#17844). In a nutshell, exome sequencing is heavily biased by the capture process and different probe affinities. A tool like VarScan2? can use a pair of matched arrays (i.e. tumor and normal) processed at the same time to do CN calls from capture data, but with sufficiently less resolution and accuracy. The bottom line is that if you want to do CN-calling, you really need whole-genome sequencing.

- May I use readDepth to calculate mouse CNV? How can I download annotation files?

Yes, you can use this package on mouse. Look at the instructions given in the Annotations section, including a link to scripts for creating your own annotations.

- What does 'absolute copy number call' mean?

Absolute copy Number call is the estimate of the number of copies at that position (as opposed to the log2 ratio). A copy-number neutral region in a diploid organism (like human) will have a value of 2.

- What does number of bins mean?

In the segmented output, the number of bins is the number of consecutive bins that were merged into a given segment.

- Can readDepth report the exact breakpoints?What are the columns of star and stop.Can I think they are positions of breakpoints?

Since this is a windowed depth based approached, the breakpoints will not be precise. CN data can be integrated with exact breakpoint calls from software such as breakdancer, using the rd.matchBreakpoints() function.

- Hi, I want to try readDepth but I only have bam files from the mapping. Does readDepth accepts bam files? If not, how can I convert .bam to .bed? Thanks

samtools view -F 4 myfile.bam | awk 'OFS="\t"{print $3,$4-1,$4}' >myfile.bed
