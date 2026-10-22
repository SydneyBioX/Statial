# Test whether an object is a kontextualResult

Test whether an object is a kontextualResult

## Usage

``` r
isKontextual(kontextualResult)
```

## Arguments

- kontextualResult:

  a object to test

## Examples

``` r
data <- data.frame()
if (!isKontextual(data)) print("Not a kontextualResult")
#> [1] "Not a kontextualResult"
```
