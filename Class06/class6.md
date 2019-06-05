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
    ## n0/dcnvrwqx7y1bmnx28rlv12q80000gp/T//RtmpAznbTy/4AKE.pdb exists. Skipping
    ## download

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## n0/dcnvrwqx7y1bmnx28rlv12q80000gp/T//RtmpAznbTy/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## n0/dcnvrwqx7y1bmnx28rlv12q80000gp/T//RtmpAznbTy/1E4Y.pdb exists. Skipping
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
rbind(s1.b, s2.b, s3.b)
```

    ##       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
    ## s1.b 29.02 18.44 16.20 19.67 20.26 20.55 17.05 22.13 26.71 33.05 30.66
    ## s2.b 37.14 25.76 23.90 17.83 19.86 21.75 20.21 16.92 17.47 18.35 18.31
    ## s3.b 25.46 17.86 10.28  4.73  4.36  5.10  9.59 12.19 11.41  9.39  8.08
    ##      [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
    ## s1.b 32.73 25.61 33.19 41.03 24.09 16.18 19.14 29.19 14.79 19.63 28.54
    ## s2.b 20.57 14.56 17.87 11.87 24.63 21.29 35.13 29.68 23.96 32.34 35.34
    ## s3.b  9.01 11.77 12.15 12.72  9.62 12.18 19.95 19.59 15.73 22.51 25.87
    ##      [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
    ## s1.b 27.49 32.56 17.13 15.50  6.98 24.07 24.00 23.94 30.70 24.70 32.84
    ## s2.b 35.64 38.91 29.00 36.55 28.83 27.15 30.28 28.13 19.90 21.95 25.07
    ## s3.b 23.08 20.97 17.28 12.69 12.24 14.14 14.05  9.38  5.03  7.78 10.13
    ##      [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44]
    ## s1.b 34.60 33.01 44.60 50.74 57.32 47.04 67.13 81.04 75.20 59.68 55.63
    ## s2.b 16.15 18.35 21.19 27.13 28.55 21.10 38.88 33.63 29.51 29.21 33.01
    ## s3.b  8.96  7.50  5.48  2.97  2.73  3.23  7.81 10.40 10.67 12.79 17.90
    ##      [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
    ## s1.b 45.12 39.04 44.31 38.21 43.70 44.19 47.00 48.67 41.54 50.22 45.07
    ## s2.b 20.92 17.17 25.84 29.80 16.89 24.66 35.62 23.52 23.37 34.41 25.96
    ## s3.b 13.56 12.94 14.78 11.31  8.79 14.13 15.10  7.92  8.15 14.28 14.04
    ##      [,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
    ## s1.b 49.77 52.04 44.82 39.75 35.79 38.92 37.93 27.18 26.86 27.53 31.16
    ## s2.b 16.79 20.20 23.72 23.29 25.23 19.81 19.00 20.21 22.62 21.40 23.47
    ## s3.b 12.42 11.84  6.57  9.59 11.84  9.61 12.18  7.89  5.74  5.31  7.67
    ##      [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77]
    ## s1.b 27.08 23.03 28.12 24.78 24.22 18.69 40.67 38.08 55.26 46.29 26.25
    ## s2.b 23.20 20.21 25.90 30.58 28.25 37.60 44.66 54.46 91.10 92.02 86.85
    ## s3.b  7.99  8.24 12.34 20.98 17.93 16.30 16.94 22.19 22.36 18.96 17.18
    ##      [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88]
    ## s1.b 37.14 27.50 16.86 27.76 19.27 22.22 26.70 25.52 21.22  15.9 15.84
    ## s2.b 80.21 68.72 42.01 27.69 23.06 21.98 18.60 20.17 15.06  14.2 23.07
    ## s3.b 18.99 16.65 13.39 11.61 10.10 11.03 13.31 12.66  9.44   6.6  5.20
    ##      [,89] [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99]
    ## s1.b 22.44 19.61 21.23 21.79 17.64 22.19 22.73 16.80 23.25 35.95 24.42
    ## s2.b 20.36 25.76 17.02 13.71 23.88 26.72 22.58 24.51 45.23 38.07 36.97
    ## s3.b  5.06  6.16  6.20  6.24  6.34  7.39  7.86 11.66 17.87 17.67 14.63
    ##      [,100] [,101] [,102] [,103] [,104] [,105] [,106] [,107] [,108] [,109]
    ## s1.b  20.96  20.00  25.99  24.39  17.19  12.16  17.35  24.97  14.08  22.01
    ## s2.b  35.17  37.83  43.69  29.14  24.56  25.20  19.27  20.88  18.27  16.96
    ## s3.b  14.30  16.98  19.84  13.36  10.93  11.52   7.56   8.85   7.07  10.08
    ##      [,110] [,111] [,112] [,113] [,114] [,115] [,116] [,117] [,118] [,119]
    ## s1.b  22.26  22.78  27.47  30.49  32.02  20.90  27.03  23.84  44.37  42.47
    ## s2.b  21.38  18.33  23.18  21.15  21.97  22.63   9.74  16.71  26.18  30.39
    ## s3.b  12.34  12.05  13.10  18.63  21.34  15.73  13.16  14.04  18.13  13.59
    ##      [,120] [,121] [,122] [,123] [,124] [,125] [,126] [,127] [,128] [,129]
    ## s1.b  33.48  44.56  56.67  60.18  66.62  59.95  70.81  88.63 100.11  86.60
    ## s2.b  22.95  25.51  20.28  16.86  21.94  20.59  21.64  27.42  35.72  23.47
    ## s3.b  12.12  13.37  10.57   6.60   7.73   7.91  11.31  14.38  14.60  12.25
    ##      [,130] [,131] [,132] [,133] [,134] [,135] [,136] [,137] [,138] [,139]
    ## s1.b  85.80  77.48  68.13  52.66  45.34  52.43  60.90  62.64  72.19  66.75
    ## s2.b  31.57  23.71  19.01  21.52  19.40  24.32  34.28  23.96  23.14  26.60
    ## s3.b  12.33  11.10  11.53  10.44   9.18  11.36  17.28  16.45  15.21  12.11
    ##      [,140] [,141] [,142] [,143] [,144] [,145] [,146] [,147] [,148] [,149]
    ## s1.b  58.73  74.57  79.29  79.53  76.58  66.40  64.76  70.48  74.84  70.11
    ## s2.b  24.94  28.49  28.18  41.64  23.85  28.67  28.76  35.16  35.46  28.74
    ## s3.b  12.12  14.10  14.94  21.72  16.82  12.61  13.40  12.64  12.24   9.13
    ##      [,150] [,151] [,152] [,153] [,154] [,155] [,156] [,157] [,158] [,159]
    ## s1.b  74.82  78.61  78.24  66.70  66.10  67.01  72.28  80.64  68.54  43.23
    ## s2.b  26.99  31.74  40.41  33.73  25.57  29.13  29.74  36.32  22.58  22.82
    ## s3.b  12.31  19.68  19.83  15.34  15.61  14.07  13.64  16.87  11.89  12.92
    ##      [,160] [,161] [,162] [,163] [,164] [,165] [,166] [,167] [,168] [,169]
    ## s1.b  51.24  45.72  61.60  45.61  42.57  41.03  41.02  33.34  19.48  34.38
    ## s2.b  46.67  29.44  25.40  17.27  20.38  21.55  19.19  15.89  18.37  30.51
    ## s3.b  19.93  23.72  23.13  13.35  11.51  18.51  17.24  11.92  12.36  13.42
    ##      [,170] [,171] [,172] [,173] [,174] [,175] [,176] [,177] [,178] [,179]
    ## s1.b  33.11  25.48  29.68  40.71  32.91  24.41  19.20  15.43  19.93  20.66
    ## s2.b  18.47  11.70  18.45  24.75  16.63  20.80  19.62  22.56  19.87  20.22
    ## s3.b  11.45  11.09  14.19  14.22  12.15  10.49  11.29  11.74   9.53   7.65
    ##      [,180] [,181] [,182] [,183] [,184] [,185] [,186] [,187] [,188] [,189]
    ## s1.b  12.72  21.40  18.21  26.68  34.50  25.77  26.52  36.85  31.05  39.84
    ## s2.b  21.16  22.13  20.66  22.82  32.86  26.04  20.60  44.44  35.28  38.03
    ## s3.b   7.21   7.56   8.14  11.07  16.93  11.12   8.79  16.03  18.87  17.72
    ##      [,190] [,191] [,192] [,193] [,194] [,195] [,196] [,197] [,198] [,199]
    ## s1.b  48.03  23.04  29.57  23.00  23.80  26.59  25.49  23.25  19.89  32.37
    ## s2.b  28.46  29.10  30.19  26.17  22.71  23.39  23.44  16.27  21.26  24.67
    ## s3.b  14.72  14.08  14.21   9.99   6.63  10.11  12.64  15.06  14.21  14.20
    ##      [,200] [,201] [,202] [,203] [,204] [,205] [,206] [,207] [,208] [,209]
    ## s1.b  30.97  42.16  29.64  29.69  33.15  26.38  23.17  29.35  32.80  25.92
    ## s2.b  19.12  23.26  21.75  24.59  27.26  22.63  26.40  31.60  29.57  30.90
    ## s3.b  16.39  16.31  16.07  17.83  20.24  14.28  17.10  17.00  18.88  17.13
    ##      [,210] [,211] [,212] [,213] [,214]
    ## s1.b  38.01  45.95  44.26  44.35  70.26
    ## s2.b  32.29  46.86  41.73  49.31  66.76
    ## s3.b  23.68  24.72  19.74  24.12  33.57

``` r
cbind(s1.b, s2.b, s3.b)
```

    ##          s1.b  s2.b  s3.b
    ##   [1,]  29.02 37.14 25.46
    ##   [2,]  18.44 25.76 17.86
    ##   [3,]  16.20 23.90 10.28
    ##   [4,]  19.67 17.83  4.73
    ##   [5,]  20.26 19.86  4.36
    ##   [6,]  20.55 21.75  5.10
    ##   [7,]  17.05 20.21  9.59
    ##   [8,]  22.13 16.92 12.19
    ##   [9,]  26.71 17.47 11.41
    ##  [10,]  33.05 18.35  9.39
    ##  [11,]  30.66 18.31  8.08
    ##  [12,]  32.73 20.57  9.01
    ##  [13,]  25.61 14.56 11.77
    ##  [14,]  33.19 17.87 12.15
    ##  [15,]  41.03 11.87 12.72
    ##  [16,]  24.09 24.63  9.62
    ##  [17,]  16.18 21.29 12.18
    ##  [18,]  19.14 35.13 19.95
    ##  [19,]  29.19 29.68 19.59
    ##  [20,]  14.79 23.96 15.73
    ##  [21,]  19.63 32.34 22.51
    ##  [22,]  28.54 35.34 25.87
    ##  [23,]  27.49 35.64 23.08
    ##  [24,]  32.56 38.91 20.97
    ##  [25,]  17.13 29.00 17.28
    ##  [26,]  15.50 36.55 12.69
    ##  [27,]   6.98 28.83 12.24
    ##  [28,]  24.07 27.15 14.14
    ##  [29,]  24.00 30.28 14.05
    ##  [30,]  23.94 28.13  9.38
    ##  [31,]  30.70 19.90  5.03
    ##  [32,]  24.70 21.95  7.78
    ##  [33,]  32.84 25.07 10.13
    ##  [34,]  34.60 16.15  8.96
    ##  [35,]  33.01 18.35  7.50
    ##  [36,]  44.60 21.19  5.48
    ##  [37,]  50.74 27.13  2.97
    ##  [38,]  57.32 28.55  2.73
    ##  [39,]  47.04 21.10  3.23
    ##  [40,]  67.13 38.88  7.81
    ##  [41,]  81.04 33.63 10.40
    ##  [42,]  75.20 29.51 10.67
    ##  [43,]  59.68 29.21 12.79
    ##  [44,]  55.63 33.01 17.90
    ##  [45,]  45.12 20.92 13.56
    ##  [46,]  39.04 17.17 12.94
    ##  [47,]  44.31 25.84 14.78
    ##  [48,]  38.21 29.80 11.31
    ##  [49,]  43.70 16.89  8.79
    ##  [50,]  44.19 24.66 14.13
    ##  [51,]  47.00 35.62 15.10
    ##  [52,]  48.67 23.52  7.92
    ##  [53,]  41.54 23.37  8.15
    ##  [54,]  50.22 34.41 14.28
    ##  [55,]  45.07 25.96 14.04
    ##  [56,]  49.77 16.79 12.42
    ##  [57,]  52.04 20.20 11.84
    ##  [58,]  44.82 23.72  6.57
    ##  [59,]  39.75 23.29  9.59
    ##  [60,]  35.79 25.23 11.84
    ##  [61,]  38.92 19.81  9.61
    ##  [62,]  37.93 19.00 12.18
    ##  [63,]  27.18 20.21  7.89
    ##  [64,]  26.86 22.62  5.74
    ##  [65,]  27.53 21.40  5.31
    ##  [66,]  31.16 23.47  7.67
    ##  [67,]  27.08 23.20  7.99
    ##  [68,]  23.03 20.21  8.24
    ##  [69,]  28.12 25.90 12.34
    ##  [70,]  24.78 30.58 20.98
    ##  [71,]  24.22 28.25 17.93
    ##  [72,]  18.69 37.60 16.30
    ##  [73,]  40.67 44.66 16.94
    ##  [74,]  38.08 54.46 22.19
    ##  [75,]  55.26 91.10 22.36
    ##  [76,]  46.29 92.02 18.96
    ##  [77,]  26.25 86.85 17.18
    ##  [78,]  37.14 80.21 18.99
    ##  [79,]  27.50 68.72 16.65
    ##  [80,]  16.86 42.01 13.39
    ##  [81,]  27.76 27.69 11.61
    ##  [82,]  19.27 23.06 10.10
    ##  [83,]  22.22 21.98 11.03
    ##  [84,]  26.70 18.60 13.31
    ##  [85,]  25.52 20.17 12.66
    ##  [86,]  21.22 15.06  9.44
    ##  [87,]  15.90 14.20  6.60
    ##  [88,]  15.84 23.07  5.20
    ##  [89,]  22.44 20.36  5.06
    ##  [90,]  19.61 25.76  6.16
    ##  [91,]  21.23 17.02  6.20
    ##  [92,]  21.79 13.71  6.24
    ##  [93,]  17.64 23.88  6.34
    ##  [94,]  22.19 26.72  7.39
    ##  [95,]  22.73 22.58  7.86
    ##  [96,]  16.80 24.51 11.66
    ##  [97,]  23.25 45.23 17.87
    ##  [98,]  35.95 38.07 17.67
    ##  [99,]  24.42 36.97 14.63
    ## [100,]  20.96 35.17 14.30
    ## [101,]  20.00 37.83 16.98
    ## [102,]  25.99 43.69 19.84
    ## [103,]  24.39 29.14 13.36
    ## [104,]  17.19 24.56 10.93
    ## [105,]  12.16 25.20 11.52
    ## [106,]  17.35 19.27  7.56
    ## [107,]  24.97 20.88  8.85
    ## [108,]  14.08 18.27  7.07
    ## [109,]  22.01 16.96 10.08
    ## [110,]  22.26 21.38 12.34
    ## [111,]  22.78 18.33 12.05
    ## [112,]  27.47 23.18 13.10
    ## [113,]  30.49 21.15 18.63
    ## [114,]  32.02 21.97 21.34
    ## [115,]  20.90 22.63 15.73
    ## [116,]  27.03  9.74 13.16
    ## [117,]  23.84 16.71 14.04
    ## [118,]  44.37 26.18 18.13
    ## [119,]  42.47 30.39 13.59
    ## [120,]  33.48 22.95 12.12
    ## [121,]  44.56 25.51 13.37
    ## [122,]  56.67 20.28 10.57
    ## [123,]  60.18 16.86  6.60
    ## [124,]  66.62 21.94  7.73
    ## [125,]  59.95 20.59  7.91
    ## [126,]  70.81 21.64 11.31
    ## [127,]  88.63 27.42 14.38
    ## [128,] 100.11 35.72 14.60
    ## [129,]  86.60 23.47 12.25
    ## [130,]  85.80 31.57 12.33
    ## [131,]  77.48 23.71 11.10
    ## [132,]  68.13 19.01 11.53
    ## [133,]  52.66 21.52 10.44
    ## [134,]  45.34 19.40  9.18
    ## [135,]  52.43 24.32 11.36
    ## [136,]  60.90 34.28 17.28
    ## [137,]  62.64 23.96 16.45
    ## [138,]  72.19 23.14 15.21
    ## [139,]  66.75 26.60 12.11
    ## [140,]  58.73 24.94 12.12
    ## [141,]  74.57 28.49 14.10
    ## [142,]  79.29 28.18 14.94
    ## [143,]  79.53 41.64 21.72
    ## [144,]  76.58 23.85 16.82
    ## [145,]  66.40 28.67 12.61
    ## [146,]  64.76 28.76 13.40
    ## [147,]  70.48 35.16 12.64
    ## [148,]  74.84 35.46 12.24
    ## [149,]  70.11 28.74  9.13
    ## [150,]  74.82 26.99 12.31
    ## [151,]  78.61 31.74 19.68
    ## [152,]  78.24 40.41 19.83
    ## [153,]  66.70 33.73 15.34
    ## [154,]  66.10 25.57 15.61
    ## [155,]  67.01 29.13 14.07
    ## [156,]  72.28 29.74 13.64
    ## [157,]  80.64 36.32 16.87
    ## [158,]  68.54 22.58 11.89
    ## [159,]  43.23 22.82 12.92
    ## [160,]  51.24 46.67 19.93
    ## [161,]  45.72 29.44 23.72
    ## [162,]  61.60 25.40 23.13
    ## [163,]  45.61 17.27 13.35
    ## [164,]  42.57 20.38 11.51
    ## [165,]  41.03 21.55 18.51
    ## [166,]  41.02 19.19 17.24
    ## [167,]  33.34 15.89 11.92
    ## [168,]  19.48 18.37 12.36
    ## [169,]  34.38 30.51 13.42
    ## [170,]  33.11 18.47 11.45
    ## [171,]  25.48 11.70 11.09
    ## [172,]  29.68 18.45 14.19
    ## [173,]  40.71 24.75 14.22
    ## [174,]  32.91 16.63 12.15
    ## [175,]  24.41 20.80 10.49
    ## [176,]  19.20 19.62 11.29
    ## [177,]  15.43 22.56 11.74
    ## [178,]  19.93 19.87  9.53
    ## [179,]  20.66 20.22  7.65
    ## [180,]  12.72 21.16  7.21
    ## [181,]  21.40 22.13  7.56
    ## [182,]  18.21 20.66  8.14
    ## [183,]  26.68 22.82 11.07
    ## [184,]  34.50 32.86 16.93
    ## [185,]  25.77 26.04 11.12
    ## [186,]  26.52 20.60  8.79
    ## [187,]  36.85 44.44 16.03
    ## [188,]  31.05 35.28 18.87
    ## [189,]  39.84 38.03 17.72
    ## [190,]  48.03 28.46 14.72
    ## [191,]  23.04 29.10 14.08
    ## [192,]  29.57 30.19 14.21
    ## [193,]  23.00 26.17  9.99
    ## [194,]  23.80 22.71  6.63
    ## [195,]  26.59 23.39 10.11
    ## [196,]  25.49 23.44 12.64
    ## [197,]  23.25 16.27 15.06
    ## [198,]  19.89 21.26 14.21
    ## [199,]  32.37 24.67 14.20
    ## [200,]  30.97 19.12 16.39
    ## [201,]  42.16 23.26 16.31
    ## [202,]  29.64 21.75 16.07
    ## [203,]  29.69 24.59 17.83
    ## [204,]  33.15 27.26 20.24
    ## [205,]  26.38 22.63 14.28
    ## [206,]  23.17 26.40 17.10
    ## [207,]  29.35 31.60 17.00
    ## [208,]  32.80 29.57 18.88
    ## [209,]  25.92 30.90 17.13
    ## [210,]  38.01 32.29 23.68
    ## [211,]  45.95 46.86 24.72
    ## [212,]  44.26 41.73 19.74
    ## [213,]  44.35 49.31 24.12
    ## [214,]  70.26 66.76 33.57

``` r
hc <- hclust( dist( rbind(s1.b, s2.b, s3.b) ) )
plot(hc)
```

![](class6_files/figure-markdown_github/unnamed-chunk-32-1.png)