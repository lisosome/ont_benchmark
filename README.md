# ONT_Benchmark

Little nexflow workflow to perform SV variant calling from BAM files of HG002 samples obtained from different alignment alghorithms and from downsampled fastqs.
The depths of coverage that will be used are:
* 5X
* 10X
* 20X
* 30X
* 40X

The alignment alghorithms used are:
* [minimap2](https://github.com/lh3/minimap2)
* [winnowmap2](https://github.com/marbl/Winnowmap)
* [lra](https://github.com/ChaissonLab/LRA)
* [ngmlr](https://github.com/philres/ngmlr)
* [HQAlign](https://github.com/joshidhaivat/HQAlign)

The SV calling algorithms are:
* [sniffles2](https://github.com/fritzsedlazeck/Sniffles)
* [CAMPHOR](https://github.com/afujimoto/CAMPHOR)
* [cuteSV](https://github.com/tjiangHIT/cuteSV)
* [Nanovar](https://github.com/cytham/nanovar)
* [Duet](https://github.com/yekaizhou/duet)
* [NanoSV](https://github.com/mroosmalen/nanosv)