Class 6: R functions
================
Kellie Lemoine
April 19, 2019

Overview
--------

Today we will focus on **R functions** but we will start with a bit of **file reading**

``` r
plot(1:10, type = "l", col="blue")
```

![](class6_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
read.table("test1.txt", header = TRUE, sep = ",")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table("test2.txt", header = TRUE, sep = "$")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table("test3.txt", header = FALSE)
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

``` r
read.csv("https://bioboot.github.io/bggn213_S19/class-material/test2.txt", header = TRUE, sep = "$")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

Our first function
------------------

Add some numbers

``` r
add <- function(x, y=1) {
  # the body
  x + y
}
```

``` r
add(4)
```

    ## [1] 5

``` r
add(4, 5)
```

    ## [1] 9

``` r
add( c(1,3,5), 1)
```

    ## [1] 2 4 6

``` r
#add( 1, 3, 5)
```

``` r
#add(x=1, y="barry")
```

our 2nd example function
------------------------

``` r
rescale <- function(x) {
 rng <-range(x)
 (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale( c(1, 3, NA, 5, 10))
```

    ## [1] NA NA NA NA NA

``` r
x <- c(1, 3, NA, 5, 10)
rng <- range(x, na.rm = TRUE)
rng
```

    ## [1]  1 10

``` r
 (x - rng[1]) / (rng[2] - rng[1])
```

    ## [1] 0.0000000 0.2222222        NA 0.4444444 1.0000000

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale2 <- function(x, na.rm = TRUE) {
 rng <-range(x, na.rm = na.rm)
 (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale2( c(1, 3, NA, 10) )
```

    ## [1] 0.0000000 0.2222222        NA 1.0000000

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {

 rng <-range(x, na.rm=TRUE)
 print("Hello")
 
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 
 print("is it me you are looking for?")

 if(plot) {
   plot(answer, typ="b", lwd=4)
   print("Dont sing please")}
 
 print("I can see it in ...")
 return(answer)
}
```

``` r
rescale3(1:10)
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"
    ## [1] "I can see it in ..."

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale3( 1:10, plot = TRUE )
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"

![](class6_files/figure-markdown_github/unnamed-chunk-20-1.png)

    ## [1] "Dont sing please"
    ## [1] "I can see it in ..."

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

section 1B of hands on worksheet
--------------------------------

``` r
#install.packages("bio3d")
```

``` r
library(bio3d)
```

``` r
s1 <- read.pdb("4AKE")
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE")
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y")
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.chainA
```

    ## 
    ##  Call:  trim.pdb(pdb = s1, chain = "A", elety = "CA")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 214,  XYZs#: 642  Chains#: 1  (values: A)
    ## 
    ##      Protein Atoms#: 214  (residues/Calpha atoms#: 214)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
    ##       DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI
    ##       VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
    ##       YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
s1.b <- s1.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor") 
```

![](class6_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
plotb3(s1.b, typ="l", ylab="Bfactor") 
```

![](class6_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## n0/dcnvrwqx7y1bmnx28rlv12q80000gp/T//RtmpFjaSkA/4AKE.pdb exists. Skipping
    ## download

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## n0/dcnvrwqx7y1bmnx28rlv12q80000gp/T//RtmpFjaSkA/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## n0/dcnvrwqx7y1bmnx28rlv12q80000gp/T//RtmpFjaSkA/1E4Y.pdb exists. Skipping
    ## download

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class6_files/figure-markdown_github/unnamed-chunk-29-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class6_files/figure-markdown_github/unnamed-chunk-29-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class6_files/figure-markdown_github/unnamed-chunk-29-3.png)

``` r
#rbind(s1.b, s2.b, s3.b)
```

``` r
#cbind(s1.b, s2.b, s3.b)
```

``` r
hc <- hclust( dist( rbind(s1.b, s2.b, s3.b) ) )
plot(hc)
```

![](class6_files/figure-markdown_github/unnamed-chunk-32-1.png)
