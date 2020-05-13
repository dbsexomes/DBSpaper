# DNA damage related metrics 


## Dependencies

* python2.7
* pysam 


## Reads with nucleotide mismatches with respect to the reference genome

```
usage: qc_mismatch.py [-h] [-i BAM] [-r BEDFILE] [-s SAMPLENAME] [-v VCFFILE] [-o OUTFILE]

	BAM	Sample alignment file
	VCFFILE	VCF file (filtered) with sample calls
	SAMPLENAME Sample's name in VCF and alignment file
	BEDFILE	Exome capture region in bed format
	OUTFILE	Output file name

 e.g. python2.7 metrics/qc_mismatch.py -i testdata/HG01606_HDAC2.bam -r refdata/nimblegen_capture.bed -s HG01606 -v testdata/1000G_chr6_114257320-114292359.vcf  -o HG01606_mismatch.tsv
```
#### Output file - important columns:


|    column		       | description                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| sample                       | sample name                                                  |
| 0_mismatch                   | reads with 0 mismatches with the reference genome            |
| 1_mismatch                   | reads with 1 mismatch with the reference genome              |
| 2_mismatch                   | reads with 2 mismatches with the reference genome            |
| 3_mismatch                   | reads with 3 or more mismatches with the reference genome    |


## Nucleotide misincorporations by base change type (NMBC)

```
usage: qc_nmbc.py [-h] [-v INVCF] [-s SAMPLEID] [-q MQMIN] [-b BQMIN]
                  [-f BAMFILE] [-r BEDFILE] [-o OUTFILE]

Compute number and fraction of non-variant base changes

optional arguments:
  -h, --help            show this help message and exit
  -v INVCF, --invcf INVCF
                        Name of VCF file with called variants
  -s SAMPLEID, --sampleid SAMPLEID
                        Sample name whose variants are to be analyzed
  -q MQMIN, --mqmin MQMIN
                        Minimum mapping quality of a read to be analyzed
  -b BQMIN, --bqmin BQMIN
                        Minimum quality of a base to be considered
  -f BAMFILE, --bamfile BAMFILE
                        Name of BAM file
  -r BEDFILE, --bedfile BEDFILE
                        Name of bed file
  -o OUTFILE, --outfile OUTFILE
                        Name of file for results

 e.g. python2.7 metrics/qc_nmbc.py -b refdata/nimblegen_capture.bed -f testdata/HG01606_HDAC2.bam -v testdata/1000G_chr6_114257320-114292359.vcf -o HG01606_nmbc.tsv
```
#### Output file - important columns:

|    column                    | description                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| sample                       | sample id                                                    |
| A>C                          | number of A>C changes in all reads/total A changes           |
| A>G                          | number of A>G changes in all reads/total A changes           |
| A>T                          | number of A>T changes in all reads/total A changes           |
| ...                          | ....                                                         |

