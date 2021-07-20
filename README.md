# viRome_legacy
Legacy R code for the R package viRome DOI:10.1093/bioinformatics/btt297

Legacy DESCRIPTION file:

```sh
Package: viRome
Type: Package
Title: Analysis and visualisation of short-read sequencing data from
        virus-infection studies
Version: 0.10
Date: 2015-11-12
Author: Mick Watson
Maintainer: Mick Watson <mick.watson@roslin.ed.ac.uk>
Depends: seqinr, plyr, gsubfn, Rsamtools, reshape2
Suggests: seqLogo, motifStack
Imports: S4Vectors
Description: During infection with a variety of viruses, host cells often raise a siRNA (~21bp) or piRNA (24-29bp) response.  The viRome package takes aligned data from a bam file and allows users to visualise and analyse the data in a number of ways
License: BSD
Packaged: 2015-11-12 10:39:36 UTC; mwatson9
```

Legacy NAMESPACE file:

```sh
exportPattern("^[[:alpha:]]+")
importClassesFrom(IRanges, Ranges, RangesList, RangedData)
importClassesFrom(S4Vectors, DataFrame)
importMethodsFrom(IRanges, append, as.matrix, as.vector,
                  as.data.frame, countOverlaps, elementLengths, end, findOverlaps,
                  follow, gsub, lapply, match, narrow, order,
                  precede, queryHits, rev, shift,
                  start, "start<-", subjectHits, 
                  unique, unlist, values, "values<-", which, width)
importFrom(IRanges, CharacterList, IRanges,
           PartitioningByWidth, successiveViews)
importFrom(S4Vectors, DataFrame, isSingleString)

```

