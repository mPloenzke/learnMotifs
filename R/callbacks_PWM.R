#' Model callback function wrapper for learning PWMs along with using fixed motifs.
#'
#' Callback for the model to output plots after a user-specified number of epochs. Model must have two input layers and utilize PWMs as filters.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords callback model
#'
#' @export
combo_PWM_Motifs <- R6::R6Class("combo_PWM_Motifs",
  inherit = KerasCallback,
  public = list(
    model = NULL,
    N = NULL,
    log_dir = NULL,
    epoch = NULL,
    plot_activations = FALSE,
    plot_filters = FALSE,
    plot_crosscorrel_motifs = FALSE,
    deNovo_data = NULL,
    jaspar_data = NULL,
    test_labels = NULL,
    test_seqs = NULL,
    num_deNovo = NULL,
    motif.names = NULL,
    filter_len = NULL,
    bg = NULL,
    pseudoprob = NULL,
    initialize = function(model = NA, N = NA, log_dir,
                          plot_activations, plot_filters, plot_crosscorrel_motifs,
                          deNovo_data, jaspar_data, test_labels, test_seqs,
                          num_deNovo, motif.names, filter_len, bg=c(.25,.25,.25,.25), pseudoprob = 1e-3) {
      self$model <- model
      self$N <- N
      self$log_dir <- log_dir
      self$epoch <- 1
      self$plot_activations <- plot_activations
      self$plot_filters <- plot_filters
      self$plot_crosscorrel_motifs <- plot_crosscorrel_motifs
      self$deNovo_data <- deNovo_data
      self$jaspar_data <- jaspar_data
      self$test_labels <- test_labels
      self$test_seqs <- test_seqs
      self$num_deNovo <- num_deNovo
      self$bg <- bg
      self$pseudoprob <- pseudoprob
      idm <- table(motif.names)[table(motif.names)>1]
      for (mm in names(idm)) {
        idx <- which(motif.names==mm)
        i <- 1
        for (mi in idx) {
          motif.names[mi] <- paste(motif.names[mi],i,sep='_')
          i <- i+1
        }
      }
      self$motif.names <- motif.names
      self$filter_len <- filter_len
    },
    on_epoch_end = function(epoch, logs = list()) {
      new_wts <- self$model$get_weights()
      new_wts[1] <- restrict_to_PWM(new_wts[[1]],bg=self$bg,pseudoprob=self$pseudoprob)
      self$model$set_weights(new_wts)
      if (self$epoch %% self$N == 0) {
        epoch_dir <- make_epoch_dir(self$log_dir,self$epoch,self$model$get_weights())
        if (self$plot_activations) {
          flat_tibble <- format_activations(self$model,self$deNovo_data,self$test_labels,input_layer=1,output_layer='deNovo_pool',buffer=0,1:self$num_deNovo)
          flat_tibble_jaspar <- format_activations(self$model,self$jaspar_data,self$test_labels,input_layer=2,output_layer='jaspar_pool',self$num_deNovo,self$motif.names)
          flat_tibble <- bind_rows(flat_tibble,flat_tibble_jaspar)
          flat_tibble$Filter <- factor(flat_tibble$Filter)
          plot_filter_activations_byClass(flat_tibble,show_beta=T,show_offset=F,fl=file.path(epoch_dir,"Filter_activations.pdf"))
        }
        if (self$plot_filters) {
          deNovo_conv <- self$model$get_weights()[[1]]
          b_conv1 <- self$model$get_weights()[[2]]
          jaspar_conv <- self$model$get_weights()[[3]]
          jaspar_b_conv1 <- self$model$get_weights()[[4]]
          deNovo_conv_list <- format_filter_motifs(deNovo_conv,type='PWM',bg=self$bg)
          jaspar_conv_list <- format_filter_motifs(jaspar_conv,type='PWM',bg=self$bg)
          plot_motifs(deNovo_conv_list,ylow=0,yhigh=2,method='custom',plotheight=NULL,fl=file.path(epoch_dir,"deNovo_filter.pdf"))
          plot_motifs(jaspar_conv_list,ylow=0,yhigh=2,method='custom',plotheight=NULL,fl=file.path(epoch_dir,"jaspar_filters.pdf"))
        }
        if (self$plot_crosscorrel_motifs) {
          activations <- get_activations(model,input_layer=1,output_layer='deNovo_conv',data=self$deNovo_data)
          dna_strings <- format_activation_motifs(activations, self$test_seqs, self$filter_len, method='alipinahi', threshold=.75)
          plot_motifs(dna_strings,ylow=0,yhigh=2,method='bits',plotheight=NULL,fl=file.path(epoch_dir,'deNovo_motifs.pdf'))
        }
      }
      self$epoch <- self$epoch + 1
    }
))
#' Model callback function wrapper for learning PWMs.
#'
#' Callback for the model to output plots after a user-specified number of epochs. Model must only have a single input layer.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords callback model
#'
#' @export
deNovo_PWM_Motifs <- R6::R6Class("deNovo_PWM_Motifs",
                            inherit = KerasCallback,
                            public = list(
                              model = NULL,
                              N = NULL,
                              log_dir = NULL,
                              epoch = NULL,
                              plot_activations = FALSE,
                              plot_filters = FALSE,
                              plot_crosscorrel_motifs = FALSE,
                              deNovo_data = NULL,
                              test_labels = NULL,
                              test_seqs = NULL,
                              num_deNovo = NULL,
                              filter_len = NULL,
                              bg = NULL,
                              pseudoprob = NULL,
                              initialize = function(model = NA, N = NA, log_dir,
                                                    plot_activations, plot_filters, plot_crosscorrel_motifs,
                                                    deNovo_data, test_labels, test_seqs,
                                                    num_deNovo,filter_len, bg=c(.25,.25,.25,.25), pseudoprob=1e-3) {
                                self$model <- model
                                self$N <- N
                                self$log_dir <- log_dir
                                self$epoch <- 1
                                self$plot_activations <- plot_activations
                                self$plot_filters <- plot_filters
                                self$plot_crosscorrel_motifs <- plot_crosscorrel_motifs
                                self$deNovo_data <- deNovo_data
                                self$test_labels <- test_labels
                                self$test_seqs <- test_seqs
                                self$num_deNovo <- num_deNovo
                                self$filter_len <- filter_len
                                self$bg <- bg
                                self$pseudoprob <- pseudoprob
                              },
                              on_epoch_end = function(epoch, logs = list()) {
                                new_wts <- self$model$get_weights()
                                new_wts[1] <- restrict_to_PWM(new_wts[[1]],bg=self$bg,pseudoprob=self$pseudoprob)
                                self$model$set_weights(new_wts)
                                if (self$epoch %% self$N == 0) {
                                  epoch_dir <- make_epoch_dir(self$log_dir,self$epoch,self$model$get_weights())
                                  if (self$plot_activations) {
                                    flat_tibble <- format_activations(self$model,self$deNovo_data,self$test_labels,input_layer=NULL,output_layer='deNovo_pool',buffer=0,1:self$num_deNovo)
                                    flat_tibble$Filter <- factor(flat_tibble$Filter)
                                    plot_filter_activations_byClass(flat_tibble,show_beta=T,show_offset=F,fl=file.path(epoch_dir,"Filter_activations.pdf"))
                                  }
                                  if (self$plot_filters) {
                                    deNovo_conv <- self$model$get_weights()[[1]]
                                    b_conv1 <- self$model$get_weights()[[2]]
                                    deNovo_conv_list <- format_filter_motifs(deNovo_conv,type='PWM',bg=self$bg)
                                    plot_motifs(deNovo_conv_list,ylow=0,yhigh=2,method='custom',plotheight=NULL,fl=file.path(epoch_dir,"deNovo_filter.pdf"))
                                  }
                                  if (self$plot_crosscorrel_motifs) {
                                    activations <- get_activations(model,input_layer=NULL,output_layer='deNovo_conv',data=self$deNovo_data)
                                    dna_strings <- format_activation_motifs(activations, self$test_seqs, self$filter_len, method='alipinahi', threshold=.75)
                                    plot_motifs(dna_strings,ylow=0,yhigh=2,method='bits',plotheight=NULL,fl=file.path(epoch_dir,'deNovo_motifs.pdf'))
                                  }
                                }
                                self$epoch <- self$epoch + 1
                              }
                            ))
