#' @title CPI_calculator.
#' @description Calculate the CPI and POP between individuals.
#' @details Calculate the CPI(combined paternity index) and POP(probability of
#'   paternity) between every two individuals in a matrix or two individuals
#'   from two matrices.
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
#' @param sexcolumn A numeric value representing the gender column, default is
#'   \code{NULL}.
#' @param threshold A numeric value representing the threshold of probability of
#'   paternity, default is 0.9999
#' @return A list contains CPI matrix, POP matrix, and parent-child duos that
#'   identified as true.

CPI_calculator <- function(x, y = NULL, genename, genenamey = NULL, genestart = 2, genesel = NULL, onecolumn = FALSE, sexcolumn = NULL, threshold = 0.9999){


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
  allele_frequency <- function(genotype_x, genotype_y){
    if(!is.null(genotype_y)){genotype <- rbind(genotype_x, genotype_y)}
    else{genotype <- genotype_x}
    genotype <- as.data.frame(genotype)

    allele_freq <- lapply(genotype, function(site){
      prop.table(table(unlist(strsplit(site, split = " "))))
    })
    return(allele_freq)
  }


  # 函数：计算累积亲权指数CPI————combined paternity index, 和亲权概率POP————probability of paternity
  paternity_index <- function(genotype_x, genotype_y, allele_freq, sex, threshold){
    # 如果y矩阵不存在，则为两个x矩阵之间计算
    if(is.null(genotype_y)){genotype_y <- genotype_x}

    # 取x矩阵的个体
    x_fun <- function(X, genotype_y, allele_freq){
      # 取y矩阵的个体
      y_fun <- function(Y, X, allele_freq){
        # 计算PI值
        PI_fun <- function(X, Y, allele){

          # 缺失值直接排除计算，PI为1
          if(is.na(X) | is.na(Y)){PI <- 1}
          else{
            X_split <- unlist(strsplit(X, split = " "))
            Y_split <- unlist(strsplit(Y, split = " "))

            # PP/PP PI = 1/p
            if(sum(X_split == Y_split) + sum(X_split == rev(Y_split)) == 4){PI <- 1/allele[intersect(X_split, Y_split)]}
            # PQ/PR PI = 1/4p
            else if(sum(X_split == Y_split) + sum(X_split == rev(Y_split)) == 1){PI <- 1/(4*allele[intersect(X_split, Y_split)])}
            # PP/PQ PI = 1/2p
            else if((sum(X_split == Y_split) == 1) & sum(X_split == rev(Y_split)) == 1){PI <- 1/(2*allele[intersect(X_split, Y_split)])}
            # PQ/PQ PI = (p+q)/4pq
            else if((((sum(X_split == Y_split) == 2) & (sum(X_split == rev(Y_split)) == 0))) | (((sum(X_split == Y_split) == 0) & (sum(X_split == rev(Y_split)) == 2)))){PI <- 0.25/allele[intersect(X_split, Y_split)[1]] + 0.25/allele[intersect(X_split, Y_split[2])]}
            # 不匹配或发生突变
            else{PI <- 0}
          }

        }
        PI <- prod(mapply(PI_fun, X, Y, allele_freq))

      }
      CPI <- apply(genotype_y, 1, y_fun, X, allele_freq)

    }
    CPI <- apply(genotype_x, 1, x_fun, genotype_y, allele_freq)
    if(is.null(genotype_y)){diag(CPI) <- 0}
    POP <- CPI/(CPI+1)

    pre_surmise <- matrix(c("Offspring_ID", "parent_ID", "CPI", "POP"), nrow = 1)
    for(i in 1:ncol(CPI)){
      pre_surmise <- rbind(pre_surmise, cbind(rep(colnames(POP)[i], length(POP[which(POP[, i] > threshold), i])), rownames(POP)[which(POP[, i] > threshold)], CPI[which(POP[, i] > threshold), i], POP[which(POP[, i] > threshold), i]))
    }

    if(!is.null(sex)){
      pre_surmise <- cbind(pre_surmise[, 1:2], sex[match(pre_surmise[, 2], sex[, 1]), 2], pre_surmise[, 3:4])
      pre_surmise[1, 3] <- "gender"
    }

    return(list(CPI = CPI, POP = POP, parchi_duo = pre_surmise))
  }


  ### /*主程序*/


  # 挑选出sex列
  if (!is.null(sexcolumn)) {sex <- as.matrix(y[, c(1, sexcolumn)])}
  else {sex <- NULL}

  # 挑选出ID列和gene列
  if (genestart > 2) {
    x <- x[, -2:-(genestart-1)]
    if (!is.null(y)) {y <- y[, -2:-(genestart-1)]}
  }
  # 格式转化
  genotype_x <- unified_format(x, genename, genesel, onecolumn)
  if(!is.null(y)){genotype_y <- unified_format(y, genenamey, genesel, onecolumn)}
  else {genotype_y <- NULL}

  # 等位基因频率
  allele_freq <- allele_frequency(genotype_x, genotype_y)

  # 计算CPI
  PTresults <- paternity_index(genotype_x, genotype_y, allele_freq, sex, threshold)

  return(PTresults)
}
