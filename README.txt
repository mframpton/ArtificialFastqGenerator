Notes for ArtificialFastqGenerator v1.0.0  --  Matthew Frampton, September 2012
-------------------------------------------------------------------------------

1. Introduction
===============

The FASTQ format [1] is the standard text-based representation for nucleotide sequences and corresponding base quality scores
that are outputted by high throughput sequencing instruments such as the Illumina Genome Analyzer. Pipelines for the analysis
of Next-Generation Sequencing (NGS) data are generally composed of a set of different publicly-available software, configured
together in order to map short reads of a genome and call variants. The fidelity of pipelines is variable, hence the ability
to evaluate their performance under different conditions would be useful. ArtificialFastqGenerator facilitates this.

ArtificialFastqGenerator takes a reference genome in FASTA format as input and outputs artificial FASTQ files. Since the
artificial FASTQs are derived from the reference, the reference provides a gold-standard for the variant calling, so enabling
evaluation of the pipeline. The user can customise DNA template/read length, gap size between paired-end reads, target
coverage, whether to use quality scores taken from existing FASTQ files, and whether to simulate sequencing errors. Detailed
coverage and error summary statistics are outputted.

ArtificialFastqGenerator is platform-independent software written in Java SE 6, and has a low memory footprint dependent on the
user-specified ``nucleobaseBufferSize'' parameter (see Section 2). On an Intel Xeon X5650 processor (12M cache, 2.66 GHz, 6.40
GT/s Intel QPI), when using real quality scores and simulating sequencing errors, and with the other parameters set to their
defaults (see Section 2), it takes about 6 hrs 20 mins to process 1 Gbase.

Any bugs should be reported to Matthew.Frampton@icr.ac.uk. If the program is used to generate data for a publication, you must
cite the following paper:

M. Frampton, R. Houlston (2012) Generation of Artificial FASTQ Files to Evaluate the Performance of Next-Generation Sequencing Pipelines.
PLoS ONE 7 (11), http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0049110 .

2. Downloading ArtificialFastqGenerator
=======================================

Unzipping ArtificialFastqGenerator.zip will give you:

- ArtificialFastqGenerator.jar (the executable jar file),
- a directory called doc containing the Javadoc,
- the license file,
- this readme file,
- files for the test described in Section 3: miniReference.fasta, test1.fastq and test2.fastq

3. Running ArtificialFastqGenerator
===================================

ArtificialFastqGenerator can be run with any reference sequence in FASTA format. Below is information about all of its user
parameters:

usage java -jar ArtificialFastqGenerator.jar -O <outputPath> -R <referenceGenomePath> -S <startSequenceIdentifier>
-F1 <fastq1ForQualityScores> -F2 <fast2ForQualityScores> -CMGCS <coverageMeanGCcontentSpread> -CMP <coverageMeanPeak>
-CMPGC <coverageMeanPeakGCcontent> -CSD <coverageSD> -E <endSequenceIdentifier> -GCC <GCcontentBasedCoverage>
-GCR <GCcontentRegionSize> -L <logRegionStats> -N <nucleobaseBufferSize> -OF <outputFormat> -RCNF <readsContainingNfilter>
-RL <readLength> -SE <simulateErrorInRead> -TLM <templateLengthMean> -TLSD <templateLengthSD> -URQS <useRealQualityScores>
-X <xStart> -Y <yStart>

-h Print usage help.
-O, <outputPath>			Path for the artificial fastq and log files, including their base name (must be specified).
-R, <referenceGenomePath> 		Reference genome sequence file, (must be specified).
-S, <startSequenceIdentifier> 		The sequence identifier in the reference after which read generation should begin (must be specified).
-F1, <fastq1ForQualityScores> 		First fastq file to use for real quality scores, (must be specified if useRealQualityScores = true).
-F2, <fast2ForQualityScores> 		Second fastq file to use for real quality scores, (must be specified if useRealQualityScores = true).
-CMGCS, <coverageMeanGCcontentSpread> 	The spread of coverage mean given GC content (default = 0.22).
-CMP, <coverageMeanPeak> 		The peak coverage mean for a region (default = 37.7).
-CMPGC, <coverageMeanPeakGCcontent>	The GC content for regions with peak coverage mean (default = 0.45).
-CSD, <coverageSD> 			The coverage standard deviation divided by the mean (default = 0.2).
-E, <endSequenceIdentifier>		The sequence identifier in the reference where read generation should stop, (default = end of file).
-GCC, <GCcontentBasedCoverage> 		Whether nucleobase coverage is biased by GC content (default = true).
-GCR, <GCcontentRegionSize>		Region size in nucleobases for which to calculate GC content, (default = 150).
-L, <logRegionStats> 			The region size as a multiple of -NBS for which summary coverage statistics are recorded (default = 2).
-N, <nucleobaseBufferSize> 		The number of reference sequence nucleobases to buffer in memory, (default = 5000).
-OF, <outputFormat>			'default': standard fastq output; 'debug_nucleobases(_nuc|read_ids)': debugging.
-RCNF, <readsContainingNfilter>		Filter out no "N-containing" reads (0), "all-N" reads (1), "at-least-1-N" reads (2), (default = 0).
-RL, <readLength> 			The length of each read, (default = 76).
-SE, <simulateErrorInRead> 		Whether to simulate error in the read based on the quality scores, (default = false).
-TLM, <templateLengthMean> 		The mean DNA template length, (default = 210).
-TLSD, <templateLengthSD> 		The standard deviation of the DNA template length, (default = 60).
-URQS, <useRealQualityScores> 		Whether to use real quality scores from existing fastq files or set all to the maximum, (default = false).
-X, <xStart> 				The first read's X coordinate, (default = 1000).
-Y, <yStart> 				The first read's Y coordinate, (default = 1000).

As stated in Section 2, ArtificialFastqGenerator.zip contains additional files which can be used to test
ArtificialFastqGenerator.jar: miniReference.fasta, test1.fastq and test2.fastq. miniReference.fasta contains about 100000
nucleobases from chromosomes 1 and 2 in the human reference genome, and 120 from chromosome 3, while test1.fastq and
test2.fastq contain 10000 paired-end reads. The command below will generate the paired-end artificial FASTQs for chromosome 1,
accepting Phred quality scores from test1.fastq and test2.fastq, and using them to simulate sequencing errors. The output path
should include the file base name (e.g. $OUTPUT_DIR/Chr1), and the -S and -E parameters are prefixes of the desired sequence
identifiers, sufficiently long to ensure a match.

java -jar ArtificialFastqGenerator.jar -R miniReference.fasta -O C1 -S ">1" -E ">" -URQS true -SE true -F1 test1.fastq -F2 test2.fastq

Apart from the artificial FASTQs, ArtificialFastqGenerator will output a file which contains the start and end indexes in the
reference sequence of all the generated reads, and a log file which contains the parameter settings and summary coverage and
error statistics (see Section 7). The user can check the log file and use FastQC to confirm that the generated FASTQs have
expected characteristics given the parameter settings.

4. Read and template DNA length
===============================

Each read in the generated FASTQs will be of the length specified by the -RL parameter. The template DNA lengths are assumed
to be normally distributed, where the -TLM parameter gives the population mean, and -TLSD, the standard deviation. Hence the
gap in number of nucleobases between two reads in a pair is: template DNA length - (2 * read length).

5. Specifying a nucleobase's coverage
=====================================

Setting a nucleobase's target coverage is highly user-customisable. To set a nucleobase's target coverage,
ArtificialFastqGenerator calculates the region's GC content, and then defines and samples from a normal distribution of
coverage levels for regions with this GC content. It calculates the distribution's mean using a Gaussian function of the GC
content. The user can customise the function by setting the coverage mean peak (the height of the bell's peak) via the CMP
parameter, the GC content at which this peak occurs (the position of the centre of the peak) via the CMPGC, and how quickly
mean coverage decreases (the width of the bell). The user can also specify the standard deviation as a multiple of the mean
via the CSD parameter, and the size of the region for which GC content is calculated via the GCR parameter.

We based the default values of these parameters on the results shown in Figure 4(b) in [2], a scatterplot of normalized
coverage of each capture probe versus its GC content for four samples sequenced by an Illumina Genome Analyzer. If they wish,
the user can obviously change the parameters in order to increase/decrease the variation in coverage across regions with
different GC content. The GCC parameter can be used to switch off the biasing of coverage based on GC content.

6. Phred quality score and error generation
===========================================

If -URQS is set to false (the default), then ArtificialFastqGenerator assigns every base in every read a high Pred score of
40, Sanger-encoding "I". The exceptions are N (unknown) reference sequence nucleobases which it gives a very low quality score
of 2, Sanger-encoding "#". If -URQS is true, then ArtificialFastqGenerator takes real quality scores from the already-existing
FASTQ files which are specified by the -F1 and -F2 parameters, (again the quality score for any N is changed to 2). It will
keep going through these FASTQs, taking the quality scores in order, until it has finished generating reads. If a generated
read is a different length to the one whose quality scores are being used, then the sequence of quality scores is
lengthened/shortened accordingly: they are lengthened by duplicating the first base quality score, and shortened by removing
the first base quality score(s). By only altering the beginning of the sequence, the trend for quality scores to deteriorate
at the end of reads is preserved.

If -SE is true, then ArtificialFastqGenerator simulates sequencing error. For each base, the error simulator calculates the
probability of an error from its Phred quality score: P_e = 10^(PHRED/-10) (see [1]). When an error is simulated, at present,
each of the 3 alternative bases have an equal probability of making the substitution. Note that errors are not simulated for
Ns.

Possible future extensions of the software are to offer the user greater choice over quality generation, and if possible, to
improve the accuracy of sequencing error by basing it on more than just the quality scores e.g. surrrounding motifs [3, 4, 5].

7. Summary coverage and error statistics
========================================

The ArtificialFastqGenerator logs regional and overall summary coverage statistics and overall error statistics. The size of
the logging regions is determined by the -LRS and -NBS parameters. An ALL nucleobase is any nucleobase in the reference
sequence, however it is marked (i.e. inc. 'N'), while an ACGT nucleobase is one marked as A/C/G/T. The regional summary
coverage statistics (1 line per region) are logged underneath a header line containing the field names. The fields are:-

(1) ALL_SO_FAR_NUM: the total number of ALL nucleobases processed so far;
(2) ALL_AV_COV: the average coverage for ALL nucleobases in this region;
(3) ALL_MIN_COV: the minimum coverage for an ALL nucleobase in this region;
(4) ALL_MAX_COV: the maximum coverage for an ALL nucleobase in this region;
(5) ACGT_NUM: the number of ACGT nucleobases in this region;
(6) ACGT_AV_COV: the average coverage for ACGT nucleobases in this region;
(7) ACGT_MIN_COV: the minimum coverage for an ACGT nucleobase in this region;
(8) ACGT_MAX_COV: the maximum coverage for an ACGT nucleobase in this region;

The overall summary coverage and error statistics are logged last underneath the header line "OVERALL STATISTICS". The overall
summary coverage statistics are equivalent to the regional statistics. The overall error statistics include the total number
of reads, nucleobase calls and ACGT nucleobase calls, and the number of these AGCT nucleobase calls for which an error was
simulated.

8. References
=============

[1] Cock P, Fields C, Goto N, Heuer M, Rice P (2009) The Sanger FASTQ file format for sequences with quality scores, and the
Solexa/Illumina FASTQ variants. Nucleic Acids Research 38:1767-71.

[2] Tewhey R, Nakano M, Wan X, Pabon-Pena C, Novak B, et al. (2009) Enrichment of sequencing targets from the human genome by
solution hybridization. Genome Biology 10: r116.

[3] Dohm J, Lottaz C, Borodina T, Himmelbauer H (2008) Substantial biases in ultra-short read data sets from high-throughput
DNA sequencing. Nucleic Acids Research 36:16 1061-1073.

[4] Nakamura K, Oshima T, Morimoto T, Ikeda S, Yoshikawa H, et al. (2011) Sequence-specific error profile of Illumina
sequencers. Nucleic Acids Research 39: e90.

[5] Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, et al. (2009) The Sequence Alignment/Map format and SAMtools.
Bioinformatics 25: 2078-2079.
