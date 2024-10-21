#' transform_data
#'
#' @param data_DTM is a DTM matrix
#'
#' @return data in the format DTM: documents as rows and words as columns
#' @export
transform_data <- function(data_DTM){

  # this function check the correct structure of the data:
  if(class(data_DTM)[1]=="DTM") data <- as.matrix(data_DTM)
  if(class(data_DTM)[1]=="data.frame"){

    if(colnames(data_DTM)[1]=="Word"){

      data <- as.data.frame(table(data_DTM$Doc, data_DTM$Word))
      data <- as.matrix(stats::reshape(data, idvar = "Var1",
                                timevar = "Var2", direction = "wide")[,-1])
    }
    if(colnames(data_DTM)[1]=="j"){

      data <- as.data.frame(table(data_DTM$i, data_DTM$j))
      data <- as.matrix(stats::reshape(data, idvar = "Var1",
                                timevar = "Var2", direction = "wide")[,-1])
    }
    if(colnames(data_DTM)[1]=="Doc"){

      data <- as.data.frame(table(data_DTM$Doc, data_DTM$Word))
      data <- as.matrix(stats::reshape(data, idvar = "Var1",
                                    timevar = "Var2", direction = "wide")[,-1])
    }
    if(colnames(data_DTM)[1]=="i"){

      data <- as.data.frame(table(data_DTM$i, data_DTM$j))
      data <- as.matrix(stats::reshape(data, idvar = "Var1",
                                timevar = "Var2", direction = "wide")[,-1])
    }

  }
  if(class(data_DTM)[1]=="matrix") stop("A DTM must be passed")#data <- as.matrix(data_DTM)

  rownames(data) <- NULL

  return(data)

}


