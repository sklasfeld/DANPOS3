Changes to script:
Samantha Klasfeld updated the script to work for python3 on 6/15/2020.
Below is the latest ReadMe file downloaded from https://sites.google.com/site/danposdoc/

Install and Usage:
* Download and extract the package danpos-#.#.#.tgz.
* Read the file "Readme.txt", or go to the the package directory danpos-#.#.# and try:
python danpos.py -h


Prerequisite:
* R version 2.9.1.
* Python 2.7, rpy2, numpy 1.5.0.
* Samtools 0.1.7 (only when input data is in sam or bam format)
* Memory ~ (genome_size/step_size) x ( replicate_count + max(2, comparison_count) ) x 8 bits.


CHANGE LOG:
Version 2.2.2:
* The function fetchValueFromWig in summits was revised to accomodate sampes that have no peak in some chromosomes.
* The function to read .Wig file now allows empty line in the file
* The Wiq function is improved

June 2013, Version 2.2.1:
*help message for Wiq and Wig2Wiq? was updated.
*the mean(), sum(), and size() functions in wig.py was updated.

April 2013, Version 2.2.0:
*The Big changes designed for DANPOS2 has finally come out with this version. DANPOS becomes a comprehensive set of tools for analyzing not only nucleosomes, but also other chromatin components such as histone modifications.

Version 2.1.4:
*minor change, some bug corrections.

Version 2.1.3:
*minor change, some bug corrections.

Feb 2013, the 2.1.2 version now:
*put back the occupancy P value, fuzziness score, and fuzziness P value in the file pooled/#.peaks.xls

Feb 2013, the 2.1.1 version now:
*comes with a dantools 0.2.0 version at http://code.google.com/p/dantools/, which could be used to do additional analysis, e.g. select peaks based on differential values or distance to promoters, or do plot around transcription start sites. Thanks to Gabriel GutiŽrrez and Jared Taylor for their suggestions on this.
*Rename the final out put file to #.allPeaks.xls and contains all nucleosomes that are either differential or not, this allow users to retrieve nucleosomes showing most or least significant changes in occupancy, fuzziness, or positions.
*provide a differential FDR value in addition to each of the previous P value, thanks to Gabriel GutiŽrrez for his suggestions on this point.

Version 2.1.0
* the path seperator '-' in previous version is now replaced by ':'
* Support .bowtie format (the default ouput format of bowtie)
* Adjust clonal signal while keep the fold change between samples.
* Give out the nucleosome fragment sizes distribution for paired-end reads too.

Version 2.0.1
* A bug in setting the step size for paired-end reads input has been corrected in the new release danpos-2.0.1.

Version 2.0.1
* A bug in setting the step size for paired-end reads input has been corrected in the new release danpos-2.0.1.
* DANPOS was used to analyze nucleosome dynamics in a paper that has just been accepted to Cell, will come online soon.

Version 2.0.0
* Support paired-end reads, fragment size correction based on pair information.
* Support input files in the very popular sam/bam format.
* Need not to specify an input format; for each input directory specified, automatively serch input files in all the supported formats (bam,sam,bed,wig).
* Come with comprehensive Comments/help information in the scripts.
* Report nucleosome gain/loss results in additional seperate files.
* Support multiple groups comparison rather than the previous pair-wise comparison.
* Support normalization of each group to a specified coverage, for the reason that some samples may have pre-known degree of global nucleosome loss relative to other samples (Gossett and Lieb, 2012).
