#' @title exclusion.
#' @description Output number of mismatched loci between every two individuals.
#' @param x A matrix or data frame containing individual names and gene
#'   information.
#' @param y Another gene information matrix or data frame, default is
#'   \code{NULL}.
#' @param genename A vector or matrix containing only all gene names.
#' @param genenamey When parameter y exists, the name of all genes in matrix y
#'   should also be specified, default is \code{NULL}.
#' @param genestart A numeric representing the first column of gene information,
#'   default is 2.
#' @param genesel A vector or matrix containing selected gene names, default is
#'   \code{NULL}.
#' @param onecolumn A logical value, if \code{TRUE}, one column represents a
#'   gene site, default is \code{FALSE}.
#' @return A matrix contains the number of mismatched loci between every two
#'   individuals.

exclusion <- function(x, y = NULL, genename, genenamey = NULL, genestart = 2, genesel = NULL, onecolumn = FALSE){

  # 函数：将基因信息转化为统一格式———— "gene1, gene2"
  unified_format <- function(gene, genename, genesel, onecolumn){
    gene <- as.matrix(gene)

    # 每个位点用一列代表的情况
    if (onecolumn) {
      genotype <- matrix(unlist(apply(gene[, -1], 2, function(s){
        lapply(strsplit(s, split = ""), function(c){
          return(paste(c[1], c[length(c)], sep = " "))
        })
      })), ncol = (ncol(gene) -1))
      rownames(genotype) <- gene[, 1]
    }
    # 每个位点用两列代表的情况
    else {
      genotype <- matrix(data = NA, nrow = nrow(gene), ncol = (ncol(gene)-1)/2)
      rownames(genotype) <- gene[, 1]
      num <- (1:(ncol(genotype)))*2
      genotype[, 1:ncol(genotype)] <- paste(gene[, num], gene[, num+1], sep = " ")
    }

    # 给每个位点加上名字
    colnames(genotype) <- genename

    # 提取出挑选的基因位点
    if(!is.null(genesel)) {
      genesel <- as.matrix(genesel)
      genotype <- genotype[, genesel]
    }

    # 将缺失值转化为NA
    genotype <- gsub("0 0", NA, genotype)

    return(genotype)
  }

  # 函数：统计不匹配个数
  mismatch <- function(genotype_x, genotype_y){
    # 如果y矩阵不存在，则为两个x矩阵之间计算
    if(is.null(genotype_y)){genotype_y <- genotype_x}

    # 取x矩阵的个体
    x_fun <- function(X, genotype_y, allele_freq){
      # 取y矩阵的个体
      y_fun <- function(Y, X, allele_freq){
        # 计算PI值
        mismatch_fun <- function(X, Y, allele){

          if(is.na(X) | is.na(Y)){return(FALSE)}
          else{
            X_split <- unlist(strsplit(X, split = " "))
            Y_split <- unlist(strsplit(Y, split = " "))

            if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 0) {return(TRUE)}
            else {return(FALSE)}
          }

        }
        mismatched_num <- sum(mapply(mismatch_fun, X, Y))

      }
      mismatched_num <- apply(genotype_y, 1, y_fun, X)

    }
    mismatched_num <- apply(genotype_x, 1, x_fun, genotype_y)
    if(is.null(genotype_y)){diag(CPI) <- NA}
    return(mismatched_num)
  }

### /*主程序*/

# 挑选出ID列和gene列
if (genestart > 2) {
  x <- x[, -2:-(genestart-1)]
  if (!is.null(y)) {y <- y[, -2:-(genestart-1)]}
}
# 格式转化
genotype_x <- unified_format(x, genename, genesel, onecolumn)
if(!is.null(y)){genotype_y <- unified_format(y, genenamey, genesel, onecolumn)}
else {genotype_y <- NULL}

# 计算不匹配个数
mismatched_num <- mismatch(genotype_x, genotype_y)

return(mismatched_num)
}
