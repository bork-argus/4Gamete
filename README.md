4Gamete analyzes DNA input to find loci pairs with all four bit patterns
between them.

```
4Gamete [-n] [-j <num>] [-g]
	-n = no output
	-j <num> = number of threads to launch -- set to number of cores in processor for optimal results.
	-g = ignore bad input lines rather than halt
```

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

Output pairs which when compared includes all four possible
bit patterns: 00, 01, 10, 11 between them.
