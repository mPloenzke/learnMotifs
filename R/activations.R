#' Calculate activation difference between classes.
#'
#' Per first-layer convolutional filter, calculate the mean activation difference between class labels.
#'
#' @param data Tibble containing columns labelled Filter, Y, Activation, Effect, and Information.
#'
#' @return Tibble with mean activation difference calculated.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords activation filter difference
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate filter ungroup distinct left_join select
#' @export
calcActivationDifference <- function(data) {
  meanactivations <- data %>%
    group_by(Y, Filter) %>%
    summarise(MeanAct=mean(Activation))
  activationDiff <- mutate(meanactivations,
         Difference=filter(meanactivations,Y==1)$MeanAct - filter(meanactivations,Y==0)$MeanAct) %>%
    ungroup() %>%
    distinct(Filter,Difference) %>%
    left_join(data %>%
                select(Filter,Effect,Information) %>%
                distinct(),
              by='Filter')
  return(activationDiff)
}

#' Calculate activations.
#'
#' Calculate first-layer convolutional filter activations.
#'
#' @param model A \code{Keras} model.
#' @param input_layer Integer value denoting the input data. If NULL, assumes a single input layer.
#' @param output_layer Integer denoting output layer or string containing the layer name.
#' @param data Data array returned from \code{\link{one_hot}} to obtain activations for.
#'
#' @return Array of activations.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords activations filter
#'
#' @importFrom keras k_function
#' @export
get_activations <- function(model,input_layer=NULL,output_layer,data) {
  if (is.null(input_layer)) {input_tensor <- model$input} else {input_tensor <- model$input[[input_layer]]}
  if (is.numeric(output_layer)) {
    activations_fn <- k_function(inputs = list(input_tensor),
                                 outputs = list(model$get_layer(index=as.integer(output_layer))$output))
  } else {
    activations_fn <- k_function(inputs = list(input_tensor),
                                 outputs = list(model$get_layer(name=output_layer)$output))
  }
  activations <- activations_fn(list(data))[[1]]
  return(activations)
}
#' Format filter activations for plotting.
#'
#' Format the first-layer filter activations for plotting with \code{\link{plot_filter_activations_byClass}}.
#'
#' @param model A \code{Keras} model.
#' @param data Data array returned from \code{\link{one_hot}} to obtain activations for.
#' @param y Vector of class labels.
#' @param input_layer Integer value denoting the input data. If NULL, assumes a single input layer.
#' @param output_layer Integer denoting output layer or string containing the layer name.
#' @param buffer Where to start extracting the beta coefficients from.
#' @param motif.names Array containing names for the filters.
#'
#' @return Tibble of formatted activations.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords activations filter format
#'
#' @importFrom tibble as.tibble
#' @importFrom dplyr bind_rows left_join mutate
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
format_activations <- function(model,data,y,input_layer=1,output_layer,buffer=0,motif.names) {
  flat1 <- as.tibble(get_activations(model,input_layer,output_layer,data[y==1,,,,drop=F]))
  flat1$y <- 1
  flat0 <- as.tibble(get_activations(model,input_layer,output_layer,data[y==0,,,,drop=F]))
  flat0$y <- 0
  flat_tibble <- bind_rows(flat1,flat0)
  if (length(motif.names) == (ncol(flat_tibble)-1)) {
    simple <- TRUE
    colnames(flat_tibble) <- c(motif.names,'Y')
  } else {
    simple <- FALSE
    num <- (ncol(flat_tibble)-1)/length(motif.names)
    colnames(flat_tibble) <- c(paste('Motif',rep(motif.names,each=num),'Loc',1:num,sep='_'),'Y')
  }
  flat_tibble$Y <- as.factor(flat_tibble$Y)
  flat_tibble <- gather(flat_tibble,Filter,Activation,1:(ncol(flat1)-1))
  flat_tibble$Filter <- factor(flat_tibble$Filter)
  W_fc1 <- as.tibble(model$get_weights()[[(length(model$get_weights())-1)]])
  b_conv1 <- as.tibble(model$get_weights()[[ifelse(is.null(input_layer),2,2*input_layer)]])
  names(W_fc1) <- 'W_fc1'
  names(b_conv1) <- 'b_conv1'
  if (simple) {
    W_fc1 <- W_fc1[(buffer+1):(buffer+length(motif.names)),1]
    W_fc1$Filter <- as.factor(motif.names)
    b_conv1$Filter <- as.factor(motif.names)
    flat_tibble <- flat_tibble %>%
      left_join(W_fc1,by='Filter') %>%
      left_join(b_conv1,by='Filter') %>%
      mutate(Dense_Activation = W_fc1*Activation)
  } else {
    W_fc1 <- W_fc1[(buffer+1):(buffer+(ncol(flat1)-1)),1]
    W_fc1$Filter <- as.factor(paste(paste('Motif',rep(motif.names,times=num),sep='_'),
                                    rep(paste('Loc',1:num,sep='_'),each=length(motif.names)),sep='_'))
    b_conv1$Filter2 <- as.factor(motif.names)
    flat_tibble$Filter2 <- as.factor(sapply(strsplit(as.character(flat_tibble$Filter),'_'),`[`, 2))
    flat_tibble <- flat_tibble %>%
      left_join(W_fc1,by='Filter') %>%
      left_join(b_conv1,by='Filter2') %>%
      mutate(Dense_Activation = W_fc1*Activation,
             Location=sapply(strsplit(as.character(flat_tibble$Filter),'_'),`[`, 4),
             Filter=Filter2) %>%
      select(-Filter2)
  }
  flat_tibble$Filter <- as.character(flat_tibble$Filter)
  return(flat_tibble)
}



