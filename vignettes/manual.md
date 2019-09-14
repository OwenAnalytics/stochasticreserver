---
title: "Stochastic Reserving Manual"
author: "R. Mark Sharp"
date: "October 03, 2017"
output: 
   - rmarkdown::html_vignette
   - rmarkdown::pdf_document
   - rmarkdown::latex_document
   - rmarkdown::word_document
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Stochastic Reserving Manual}
  %\usepackage[UTF-8]{inputenc}
---


[Introduction]  
[Summary of Major Functions]   
[Installation]  
[Input]    
[Summary Statistics]    
[Software Development]    


## Introduction  
Maximum likelihood estimators provide a powerful statistical tool. In this paper we directly deal with non-linear reserving models, without the need to transform those models to make them tractable for linear or generalized linear methods. We also show how the same general approach can be easily adapted to provide estimates for a very wide range of reserving methods and models, making use of the same framework, and even much of the same computer code. We focus on the triangle of incre- mental average costs, and show how five common methods can be set in a stochastic framework.

**For more information see:**
    [A Flexible Framework for Stochastic Reserving Models
    ](http://www.variancejournal.org/issues/07-02/123.pdf "Original Paper")

## Summary of Major Functions  

### Quality Control
stub

## Installation

You can install **stochasticreserver** from github with:


```r
install.packages("devtools")
devtools::install_github("rmsharp/stochasticreserver")
```

All missing dependencies should be automatically installed.



## Input  
stub


## Summary Statistics  
stub


## Software Development

stub

