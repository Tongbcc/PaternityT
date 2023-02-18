#' @title Forensic_calculator.
#' @description Calculate the forensic parameters for every site in the group.
#' @details Forensic parameters including MAF:minor allele frequency, H_obs:
#'   observed heterozygosity, H_exp: expected heterozygosity, Ae: effective
#'   number of allele, DP: discrimination power, PE: probability of exclusion,
#'   PIC:polymorphism information content
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
#' @return A matrix contains site names and 7 forensic parameters.

Forensic_calculator <- function(x, y = NULL, genename, genenamey = NULL, genestart = 2, genesel = NULL, onecolumn = FALSE){


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


  # 函数：计算等位基因频率allele_frequency和基因型频率genotype_frequency
  allele_frequency <- function(genotype_x, genotype_y){
    if(!is.null(genotype_y)){genotype <- rbind(genotype_x, genotype_y)}
    else{genotype <- genotype_x}
    genotype <- as.data.frame(genotype)

    allele_freq <- lapply(genotype, function(site){
      prop.table(table(unlist(strsplit(site, split = " "))))
    })

    genotype_freq <- lapply(genotype, function(site){
      prop.table(table(site))
    })
    return(list(allele_freq = allele_freq, genotype_freq = genotype_freq))
  }


  # 函数：计算次等位基因频率MAF
  MAF <- function(allele_freq){
    MAF <- lapply(allele_freq, function(a){
      a[order(a, decreasing = T)][2]
    })
  }


  # 函数：计算观测杂合度H_obs和期望杂合度H_exp
  H <- function(genotype_freq, allele_freq){
    H_obs <- lapply(genotype_freq, function(a){
      b <- strsplit(names(a), split = " ")
      c_fun <- function(b, a){
        if(b[1] == b[2]){return(0)}
        else{return(a)}
      }
      H_obs <- sum(mapply(c_fun, b, a))
    })

    a_fun <- function(a, num){
      (num * (1 - sum(a^2)))/(num - 1)
    }
    H_exp <- lapply(allele_freq, FUN = a_fun, length(allele_freq))

    return(list(H_obs = H_obs, H_exp = H_exp))
  }


  # 函数：计算有效等位基因数Ae
  Ae <- function(H_obs){
    lapply(H_obs, function(a){
      Ae <- 1/(1-a)
    })
  }


  # 函数：计算个体识别能力DP
  DP <- function(allele_freq){
    lapply(allele_freq, function(a){
      1 - (2*(sum(a^2)^2)) + sum(a^4)
    })
  }


  # 函数：计算非父排除概率PE
  PE <- function(allele_freq){
    PE_0 <- lapply(allele_freq, function(a){
      PE_1 <- sum(a^2 * (1-a)^2)
      PE_2 <- sum(apply(combn(a, 2), 2, function(b){
        2 * b[1] * b[2] * (1 - sum(b))^2
      }))
      PE_0 <- PE_1 + PE_2
      return(PE_0)
    })

    return(PE_0)
  }


  # 函数：计算多态信息含量PIC
  PIC <- function(allele_freq){
    PIC <- lapply(allele_freq, function(a){
      1 - sum(a^2) - sum(a^2)^2 + sum(a^4)
    })
    return(PIC)
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

  # 等位基因频率和基因型频率
  freq <- allele_frequency(genotype_x, genotype_y)
  allele_freq <- freq[["allele_freq"]]; genotype_freq <- freq[["genotype_freq"]]

  # 计算MAF
  MAF <- MAF(allele_freq)

  # 计算杂合度
  H <- H(genotype_freq, allele_freq)
  H_obs <- H[["H_obs"]]; H_exp <- H[["H_exp"]]

  # 计算有效等位基因数
  Ae <- Ae(H_obs)

  # 计算个体识别能力DP
  DP <- DP(allele_freq)

  # 计算非父排除概率PE
  PE <- PE(allele_freq)

  # 计算多态信息含量PIC
  PIC <- PIC(allele_freq)

  final_selection <- cbind(colnames(genotype_x), MAF, H_obs, H_exp, Ae, DP, PE, PIC)
  colnames(final_selection)[1] <- "site_name"

  return(final_selection)
}
