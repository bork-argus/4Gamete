4Gamete analyzes binary bi-allelic Single Nucleotide Polymorphism (SNP) input 
to find loci pairs with all four bit patterns between them.

```
4Gamete [-n] [-j <num>] [-g]
	-n = no output
	-j <num> = number of threads to launch -- set to number of cores in processor for optimal results.
	-g = ignore bad input lines rather than halt
```

Expects files with a binary SNP array in the form of loci as rows (with row 
name column) and samples as columns.  Tab or Space delim.

Input files look like:
```
BE017350.1_124 0 0 0 0 0 0 0 1 1 1 0 0 1 1 0 1 0 0 0 1 0 1 0 1
BE017350.1_248 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 1 0 1 0 1
BE017350.1_368 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 1 0 0 0 1 0 1 0 1
BE017350.1_635 1 1 1 1 1 1 0 1 1 0 0 1 0 0 1 0 1 1 1 0 1 0 0 0
BE017350.1_653 0 1 1 1 1 0 0 0 0 0 0 1 0 0 1 0 0 1 1 0 1 0 0 0
BE017350.1_674 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 0 1 1 1 0 1 0 1 0
BE017350.1_719 0 1 0 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 1 0 0 0
BE017350.1_737 0 0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
BE017350.1_767 0 1 0 0 1 0 0 0 0 0 1 1 0 0 1 0 0 1 1 0 1 0 1 0
BE017350.1_839 1 1 1 1 1 1 0 0 0 0 1 1 0 0 1 0 1 1 1 0 1 0 1 0
```
Rows with missing data (. or *) will halt the program unless the `-g` flag is 
used; in which case that locus will not be processed.

This format can easily be achieved using the program bcftools (not included):

`bcftools query -f '%CHROM\_%POS[ %GT]' input.vcf.gz > binary_alleles_output.txt`

Output pairs which when compared includes all four possible
bit patterns: 00, 01, 10, 11 between them.
