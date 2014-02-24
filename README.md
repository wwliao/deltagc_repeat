# deltagc_repeat
**deltagc_repeat** is a Python script calculating and comparing the difference 
of GC ratio between the siRNA-accumulated regions and the non-siRNA-accumulated
 regions in each tandem repeats and randomized tandem repeats.

## How it works
Within tandem repeats, the sequences are further classified into two regions, 
with or without siRNA accumulation, for the GC-content analysis. The GC-content
percentage (GC%) based on the genomic sequence is defined as:

GC% = (G + C) / (A + T + C + G) * 100%

For each repeat, we calculated the difference of GC ratio between the 
siRNA-accumulated regions (GC<sub>siRNA</sub>%) and the non-siRNA-accumulated regions 
(GC<sub>non</sub>%), as 

deltaGC = GC<sub>siRNA</sub>% - GC<sub>non</sub>% 

We compared the deltaGCs in repeats with the randomized genomic sequences, in 
which the genomic sequences of tandem repeats are shuffled.

## How to use
Please see the command-line help:
    python deltagc_repeat.py -h

## Demo time
    python deltagc_repeat.py demo/siRNA_hits_chr1.txt demo/tandem_repeats_with_siRNA_chr1.txt demo/TAIR10_chr1.fa

The output box plot:
![](https://github.com/wwliao/deltagc_repeat/blob/master/demo/deltagc_boxplot.png)

## Where to get it
The source code is currently hosted on GitHub at: 
https://github.com/wwliao/deltagc_repeat

## Dependencies
**deltagc_repeat** needs Python 2.7 and the packages listed below:
- [Biopython](www.biopython.org/): for parsing FASTA files
- [SciPy](http://www.scipy.org): for calculating the T-test scores
- [Matplotlib](http://matplotlib.sourceforge.net/): for plotting


## License
MIT License
