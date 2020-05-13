# Variant quality related metrics 


## Dependencies

* python2.7


## Confident sites 

```
usage: qc_sites.py <input gvcf file> <output> <bedfile>

	input   Input GVCF file for sample
	output  Output file name
	bedfile	Exome capture region or other region of interest in bed format

 e.g. python2.7 metrics/qc_sites.py -i testdata/HG01606_HDAC2.g.vcf HG01606_HDAC2_sites.tsv refdata/HDAC2_capture_targets.bed
```
#### Output file - important columns:

|    column		       | description                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| sample                       | sample name                                                  |
| all_sites		       | total sites							|			
| conf_ref			| confident ref sites (GQ>=30)					|
| lq_ref			| low quality ref sites (GQ<30)					|
| conf_snp			| confident variant sites (GQ>=30)				|
| lq_snp			| low quality variant sites (GQ<30				|
| no_call			| no call site							|
| hq_fr				| fraction confident ref and variant sites			|
| lq_fr                   	| fraction low confident ref and variant sites                  |


## Variant counts (from VCF) 

* High quality variants
* Likely false positive variants (LFPV)

#### Prerequisites
* Filtered VCF file
* 1000 genomes frequence annotation in INFO field of the VCF file with format KGDB_AF
* GT and GQ fields available for calls 

```
usage: qc_var.py [-h] [-i INFILE] [-b BEDFILE] [-q MINGQ] [-f FREQ]
                 [-o OUTFILE]

Variant statistics in the capture region

optional arguments:
  -h, --help  show this help message and exit
  -i INFILE   input vcf file
  -b BEDFILE  Bed file
  -q MINGQ    min genotype quality cutoff (default 30.0)
  -f FREQ     freq cutoff
  -o OUTFILE  output file

 e.g. python2.7 metrics/qc_var.py -b refdata/nimblegen_capture.bed -i testdata/kg3_HDAC2_annot.vcf -o test_var.tsv
```
#### Output file - important columns:
|    column                    | description                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| sample                       | sample id                                                    |
| snps_hq_freq			| High quality common SNVs (freq >= 0.001 and PASS, GQ>=30	|
| snps_hq_rarenovel		| High quality rare and novel SNVs (freq < 0.001 and PASS, GQ>=30|
| snps_hq_novel2		| Likely false positive high quality SNVs			|
| indels_hq_freq		| High quality common Indels (freq >= 0.001 and PASS, GQ>=30	|
| indels_hq_rarenovel		| High quality rare and novel Indels (freq < 0.001 and PASS, GQ>=30|
| indels_hq_novel2		| Likely false positive high quality Indels			|

## Transition/Transversion ratios 

#### Prerequisites
* Filtered VCF file
* 1000 genomes frequence annotation in INFO field of the VCF file with format KGDB_AF
* GT and GQ fields available for calls


usage: qc_tstv.py [-h] [-i INFILE] [-b BEDFILE] [-q MINGQ] [-f FREQ] [-d]
                  [-o OUTFILE]


optional arguments:
  -h, --help  show this help message and exit
  -i INFILE   input vcf file
  -b BEDFILE  Bed file
  -q MINGQ    min genotype quality cutoff (default 30.0)
  -f FREQ     freq cutoff
  -d          detailed report
  -o OUTFILE  output file

 e.g. python2.7 metrics/qc_tstv.py -b refdata/nimblegen_capture.bed -i testdata/kg3_HDAC2_annot.vcf -o test_tstv.tsv

#### Output file - important columns:
|    column                    | description                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| sample                       | sample id                                                    |
| tstv_all                     | Transition/Transversion for all hiqh quality SNVs (PASS, GQ>=30)      |
| tstv_fr                     | Transition/Transversion for all common hiqh quality SNVs (freq >= 0.001 and PASS, GQ>=30)      |
| tstv_rarenovel              | Transition/Transversion for all rare hiqh quality SNVs (freq < 0.001 and PASS, GQ>=30)      |



## High quality SNV Nucleotide change frequencies 

#### Prerequisites
* Filtered VCF file
* 1000 genomes frequence annotation in INFO field of the VCF file with format KGDB_AF
* GT and GQ fields available for calls

```
usage: qc_nc.py [-h] [-i INFILE] [-b BEDFILE] [-q MINGQ] [-f FREQ]
                [-o OUTFILE]


optional arguments:
  -h, --help  show this help message and exit
  -i INFILE   input vcf file
  -b BEDFILE  Bed file
  -q MINGQ    min genotype quality cutoff (default 30.0)
  -f FREQ     freq cutoff
  -o OUTFILE  outfile name

e.g. python2.7 metrics/qc_nc.py -i testdata/kg3_HDAC2_annot.vcf -b refdata/nimblegen_capture.bed -o test_nc.tsv
```
#### Output file - important columns:
|    column                    | description                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| sample                       | sample id                                                    |
| A>C				|  A>C changes/total high quality SNVs (PASS and GQ>=30) |
| A>G				|  A>G changes/total high quality SNVs (PASS and GQ>=30) |
| A>T				|  A>T changes/total high quality SNVs (PASS and GQ>=30) |
| ...	| ... |
