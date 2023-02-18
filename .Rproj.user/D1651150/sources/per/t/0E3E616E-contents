#' @title CFSI_calculator.
#' @description Calculate the CFSI and CIBS between individuals.
#' @details Calculate the CFSI(combined full sibling index) and CIBS(combined
#'   identity by state score) between every two individuals in a matrix or two
#'   individuals from two matrices.
#' @param x A matrix or data frame containing individual names and gene
#'   information.
#' @param genename A vector or matrix containing only all gene names.
#' @param genestart A numeric representing the first column of gene information,
#'   default is 2.
#' @param genesel A vector or matrix containing selected gene names, default is
#'   \code{NULL}.
#' @param onecolumn A logical value, if \code{TRUE}, one column represents a
#'   gene site, default is \code{FALSE}.
#' @return A list contains CFSI matrix and CIBS matrix.

CFSI_calculator <- function(x, genename, genestart = 2, genesel = NULL, onecolumn = FALSE){
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


  # 函数：计算等位基因频率allele_frequency
  allele_frequency <- function(genotype_x){
    genotype <- as.data.frame(genotype_x)

    allele_freq <- lapply(genotype, function(site){
      prop.table(table(unlist(strsplit(site, split = " "))))
    })
    return(allele_freq)
  }


  # 函数：计算FSI
  full_sibling_index <- function(genotype_x, allele_freq){
    # 计算同一个矩阵中每两个个体间的CFSI
    genotype_y <- genotype_x

    # 取x矩阵的个体
    x_fun <- function(X, genotype_y, allele_freq){
      # 取y矩阵的个体
      y_fun <- function(Y, X, allele_freq){
        # 计算PI值
        FSI_fun <- function(X, Y, allele){

          # 缺失值直接排除计算，PI为1
          if(is.na(X) | is.na(Y)){FSI <- 1}

          else{
            X_split <- unlist(strsplit(X, split = " "))
            Y_split <- unlist(strsplit(Y, split = " "))

            # PP/PP
            if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 4){FSI <- 0.5 + 0.5/allele[intersect(X_split, Y_split)]}
            # PP/QQ;PP/QR;PQ/RR;PQ/RS
            else if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 0){FSI <- 0.25}
            # PQ/PR;PR/QR
            else if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 1){FSI <- 0.25 + 0.125/allele[intersect(X_split, Y_split)]}
            # PP/PQ;PQ/QQ以及PQ/PQ
            else if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 2){
              # PP/PQ;PQ/QQ
              if(X_split[1] == X_split[2] | Y_split[1] == Y_split[2]){FSI <- 0.25 + 0.25/allele[intersect(X_split, Y_split)]}
              # PQ/PQ
              else if((sum(X_split == Y_split) == 2) | (sum(X_split == rev(Y_split)) == 2)){FSI <- 0.25 + (sum(allele[intersect(X_split, Y_split)]) + 1)/(8*prod(allele[intersect(X_split, Y_split)]))}
            }
          }

        }
        FSI <- prod(mapply(FSI_fun, X, Y, allele_freq))

      }
      CFSI <- apply(genotype_y, 1, y_fun, X, allele_freq)

    }
    CFSI <- apply(genotype_x, 1, x_fun, genotype_y, allele_freq)
    if(is.null(genotype_y)){diag(CPI) <- 0}

    return(CFSI)
  }

  # 函数：计算CIBS
  ibs_cal <- function(genotype_x){
    # 计算统一矩阵中每两个个体间的CIBS
    genotype_y <- genotype_x

    # 取x矩阵的个体
    x_fun <- function(X, genotype_y){
      # 取y矩阵的个体
      y_fun <- function(Y, X){
        # 计算PI值
        IBS_fun <- function(X, Y){

          # 缺失值直接排除计算，PI为1
          if(is.na(X) | is.na(Y)){IBS <- 0}

          else{
            X_split <- unlist(strsplit(X, split = " "))
            Y_split <- unlist(strsplit(Y, split = " "))

            # PP/PP
            if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 4){IBS <- 2}
            # PP/QQ;PP/QR;PQ/RR;PQ/RS
            else if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 0){IBS <- 0}
            # PQ/PR;PR/QR
            else if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 1){IBS <- 1}
            # PP/PQ;PQ/QQ以及PQ/PQ
            else if(sum(sum(X_split == Y_split), sum(X_split == rev(Y_split))) == 2){
              # PP/PQ;PQ/QQ
              if(X_split[1] == X_split[2] | Y_split[1] == Y_split[2]){IBS <- 1}
              # PQ/PQ
              else if((sum(X_split == Y_split) == 2) | (sum(X_split == rev(Y_split)) == 2)){IBS <- 2}
            }
          }

        }
        IBS <- sum(mapply(IBS_fun, X, Y))

      }
      CIBS <- apply(genotype_y, 1, y_fun, X)

    }
    CIBS <- apply(genotype_x, 1, x_fun, genotype_y)
    if(is.null(genotype_y)){diag(CPI) <- NA}

  }


  ### /*主程序*/

  # 挑选出ID列和gene列
  if (genestart > 2) {
    x <- x[, -2:-(genestart-1)]
    if (!is.null(y)) {y <- y[, -2:-(genestart-1)]}
  }
  # 格式转化
  genotype_x <- unified_format(x, genename, genesel, onecolumn)

  # 等位基因频率
  allele_freq <- allele_frequency(genotype_x)

  # 计算CFSI
  CFSI <- full_sibling_index(genotype_x, allele_freq)

  # 计算CIBS
  CIBS <- ibs_cal(genotype_x)

  return(list(CFSI = CFSI, CIBS = CIBS))
}












