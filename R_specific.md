## Formats and Styles

**Removing # from the markdown outputs**
````r

 ```{r include=FALSE}
    knitr::opts_chunk$set(comment = NA)
 ```
````
Source : [StackOverflow](https://stackoverflow.com/questions/15081212/remove-hashes-in-r-output-from-r-markdown-and-knitr)

---
**Printing tables nice and tidy**
```r
df%>%
  kable %>%
  kable_styling()
  
```
source: [kable package](https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html)
