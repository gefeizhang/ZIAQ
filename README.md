ZIAQ: A quantile regression method for differential expression analysis of single-cell RNA-seq data
===============
A statistical differential expression analysis method that accounts for both dropout rates and complex scRNA-seq data distributions in the same model

Installation Instructions
------------
You may directly download source file `ZIAQ_1.0.tar.gz` from github link  
https://github.com/gefeizhang/ZIAQ/blob/master/ZIAQ_1.0.tar.gz

Or you may use 'devtools' to install it from github

    install_github("gefeizhang/ZIAQ")

Examples
------------
Run `ZIAQ` for indiviual gene

    y = round(100* runif(100))
    colDat = data.frame(condition = rep(c(1, 0), e = 50))
    res = ziaq_fit(y, colDat = colDat,  formula = ~ condition,
              group = 'condition', probs = c(0.25, 0.5, 0.75),
              log_i = TRUE )

Run `ZIAQ` for scRNA gene matrix
    
    # simulate gene matrix
    ymatrix = matrix(round(100* runif(100*500)), ncol = 100)
    rownames(ymatrix) = paste0('gene', 1:500)
    
    # simulate cell conditions
    colDat = data.frame(condition = rep(c(1, 0), e = 50))
    
    res = ziaq(ymatrix, colDat, formula = ~ condition,
          group = 'condition', probs = c(0.25, 0.5, 0.75),
          log_i = TRUE, parallel = FALSE, no.core = 1)


Citation
----------------
[ZIAQ: A quantile regression method for differential expression analysis of single-cell RNA-seq data](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa098/5735412?guestAccessKey=38e5976d-c09a-4138-8a56-92810ea04e8d)
Wenfei Zhang, Ying Wei, Donghui Zhang, Ethan Y Xu, Bioinformatics, btaa098
