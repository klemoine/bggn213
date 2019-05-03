Class 7: R functions and packages
================
Kellie Lemoine
April 24, 2019

More on function writing
------------------------

First we will revisit out function from last class

``` r
source("http://tinyurl.com/rescale-R") 
```

Test the ""rescale()\*\* function

``` r
rescale( 1:10 )
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale( c(1:10, "string"))
```

``` r
x <-  c(1:10, "string")
!is.numeric(x)
```

``` r
##rescale2(x)
```

Function practice
-----------------

Write a function to identify NA elements in two vectors

Start with a simple example input where I know what the answer should be

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

I am looking for the positions where it is TRUE in both vectors...

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take the sum to find how many

``` r
sum( is.na(x) & is.na(y) )
```

    ## [1] 1

This is my working snippet of code that I can use as the body of my first function

``` r
both_na <- function(x, y) {
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(x, y)
```

    ## [1] 1

``` r
both_na(c(NA, NA, NA), c(NA, NA, 1))
```

    ## [1] 2

``` r
both_na(c(NA, NA, NA), c(1, NA, NA))
```

    ## [1] 2

``` r
both_na(c(NA, NA, NA), c(1, NA, NA, NA))
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
both_na(c(NA, NA, NA), c(1, NA, NA, NA, NA, NA))
```

    ## [1] 5

Check the length of our inputs are equal

``` r
x <- c(NA, NA, NA)
y <- c(1, NA, NA, NA, NA, NA)
length(x) != length(y)
```

    ## [1] TRUE

``` r
3 != 3
```

    ## [1] FALSE

Try thr both\_na3() function with extra features

``` r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

both_na3(x, y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

The which function tells you which position is true

``` r
which(c(F, F, T))
```

    ## [1] 3

Write function to grade students homework with the lowest grade dropped
-----------------------------------------------------------------------

``` r
student1 <- c(100,100, 100, 100, 100, 100, 100, 90)
```

``` r
lowest.hw <- min(student1)
lowest.hw
```

    ## [1] 90

``` r
sum.student1 <- sum(student1) - lowest.hw
sum.student1
```

    ## [1] 700

``` r
score.student1 <- sum.student1/(length(student1)-1)
score.student1
```

    ## [1] 100

Professor's method:
-------------------

``` r
x <- c(100,100, 100, 100, 100, 100, 100, 90)

grade <- function(x) {
  (sum(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / (length(x) - 1)
}
```

``` r
#grade(x)
#grade(student2)
```

Determine grades for the entire class
-------------------------------------

``` r
url <- "http://tinyurl.com/gradeinput"

students <- read.csv(url, row.names = 1)
head(students)
```

    ##           hw1 hw2 hw3 hw4 hw5
    ## student-1 100  73 100  88  79
    ## student-2  85  64  78  89  78
    ## student-3  83  69  77 100  77
    ## student-4  88  NA  73 100  76
    ## student-5  88 100  75  86  79
    ## student-6  89  78 100  89  77

``` r
grade(students[1,])
```

    ## [1] 91.75

``` r
ans <- apply(students, 1, grade)
```

``` r
sort(ans, decreasing = TRUE)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

One last function example
-------------------------

Find the intersection of two sets

``` r
x <- df1$IDs
y <- df2$IDs

intersect(x, y)
```

    ## [1] "gene2" "gene3"

``` r
x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
y %in% x
```

    ## [1]  TRUE FALSE  TRUE FALSE

``` r
gene_intersect <- function(x, y) {
  cbind( x[x %in% y],
        y[y %in% x] )
}
```

``` r
merge(df1, df2, by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

``` r
## install.packages("BiocManager")
```
