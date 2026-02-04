#' Training the ropls model
#'
#' @param X n_samples x n_features matrix
#' @param y vector of labels
#'
#' @return ropls object
#' @export
train_ropls <- function(X, y, options = list()) {
  # suppress annoying "error"s from ropls
  #sink(file = tempfile())
  if ("n_LV" %in% names(options)) {
    predI <- options$n_LV
  } else {
    predI <- NA
  }
  try_out <- try(model <- ropls::opls(X, y,
                                      #crossValI = 5, # TODO make this an option
                                      permI = 0, # no permutation and other output to save computation time
                                      predI = predI,
                                      info.txtC = "none",
                                      fig.pdfC = "none"#,
                                      #silent = TRUE
  )
  )
  if (is(try_out, "try-error")) {
    # No model could be build, provide a model with one latent variable
    # (even if this is not significant)
    model <- ropls::opls(X, y, predI = 1,
                         #crossValI = 5,
                         permI = 0,
                         info.txtC = "none",
                         fig.pdfC = "none"#,
                         #silent = TRUE
    )
  }

  # print error messages again
  #sink()
  return(model)
}

#' Prediction using a ropls object
#'
#' @param model ropls object
#' @param X n_samples x n_features matrix
#'
#' @return vector of predicted labels
#' @export
predict_ropls <- function(model, X) {
  # suppress annoying "error"s from ropls
  #sink(file = tempfile())
  y_pred <- ropls::predict(model, newdata = X)#, silent = TRUE)
  #sink()
  return(y_pred)
}

#' PCA using ropls
#'
#' @param X n_samples x n_features matrix
#' @return ropls object
#'
#' @export
pca_ropls <- function(X) {
  # suppress annoying "error"s from ropls
  #sink(file = tempfile())
  try_out <- try(
    model <- ropls::opls(X, predI = NA,
                         #crossValI = 5, # TODO make this an option
                         permI = 0, # no permutation and other output to save computation time
                         info.txtC = "none",
                         fig.pdfC = "none"#,
                         #silent = TRUE
    )
  )
  if (is(try_out, "try-error") | ropls::getSummaryDF(model)$pre < 2) {
    # to ensure that the model has at least two prinicipal components
    model <- ropls::opls(X, predI = 2,
                         #crossValI = 5,
                         permI = 0,
                         info.txtC = "none",
                         fig.pdfC = "none"#,
                         #silent = TRUE
    )
  }
  #sink()
  return(model)
}

#' Accuracy
#'
#' @param y true labels
#' @param y_pred predicted labels
#'
#' @return classification accuracy
#' @export
score_accuracy <- function(y, y_pred) {
  num <- as.numeric(y)
  num_pred <- as.numeric(y_pred)
  correct <- which(num == num_pred)
  accuracy <- length(correct) / length(y)
  return(accuracy)
}


#' Mean squared error
#'
#' @param y true labels
#' @param y_pred predicted labels
#'
#' @return mean squared error
#' @export
score_mse <- function(y, y_pred) {
  num <- as.numeric(y)
  num_pred <- as.numeric(y_pred)
  mse <- mean((num - num_pred) ^ 2)
  return(mse)
}


#' R^2 (coefficient of determination)
#'
#' @param y true labels
#' @param y_pred predicted labels
#'
#' @return R2 value
#' @export
score_r2 <- function(y, y_pred) {
  num <- as.numeric(y)
  num_pred <- as.numeric(y_pred)
  y_bar <- mean(num)
  total_var <- sum((num - y_bar) ^ 2)
  residual_var <- sum((num - num_pred) ^ 2)
  r2 <- 1 - residual_var / total_var
  return(r2)
}

#' P (Pearson correlation)
#'
#' @param y true labels
#' @param y_pred predicted labels
#'
#' @return P
#' @export
score_pearson <- function(y, y_pred) {
  tmp <- cor.test(y, y_pred, method = "pearson")
  P <- tmp$estimate
  pvals_label <- tmp$p.value
  return(P)
}

#' P (Spearman correlation)
#'
#' @param y true labels
#' @param y_pred predicted labels
#'
#' @return P
#' @export
score_spearman <- function(y, y_pred) {
  tmp <- cor.test(y, y_pred, method = "spearman")
  P <- tmp$estimate
  pvals_label <- tmp$p.value
  return(P)
}

#' Cross-validating an mPLSDA model
#'
#' @param X n_samples x n_features matrix
#' @param y vector of labels
#'
#' @return ropls object
#' @export
#'
cross_validation_unpaired <- function(X, y, method, options, n_trials) {
  X <- as.matrix(X)
  y <- y

  vals_list <- list()

  if ("select" %in% names(method) ) {
    select <- method$select
  } else {
    select <- function(X,y) {return(colnames(X))}
    if (!("rf_trials" %in% names(options))) {
      options$rf_trials <- 0
    } else if (options$rf_trials != 0) {
      message("Warning in validate():")
      message("    no feature selector given but rf_trials != 0")
      message("    forcing rf_trials = 0")
      options$rf_trials <- 0
    }
  }

  # default to five-fold cross-validation
  if (!("n_folds" %in% names(options))) {
    options$n_folds <- 5
  }

  # also give these guys some shorter names
  train <- method$train
  predict <- method$predict
  score <- method$score


  for (trial in 1:n_trials) {
    message(paste("validate_repeat: trial", trial, "/", n_trials))

    # ----------------- INITIAL PROCESSING ----------------- #

    # see if a feature selector is passed, otherwise default
    # to selecting all features. in the latter case, also make
    # sure rf_trials = 0 since random features don't make sense



    # for permuted labels, to which the predicted outcome should be compared to with the score function
    if (!("compare_pred" %in% names(options))) {
      options$compare_pred <- "y"
    } else {
      if (!(options$compare_pred %in% c("y", "y_perm"))) {
        stop("options$compare_pred needs to be \"y\" or \"y_perm\"")
      }
    }

    # for paired data, take structure into account
    if (!("paired" %in% names(options))) {
      options$paired <- FALSE
    }
    if (options$paired) {
      if (!("X_label" %in% names(options))) {
        stop("X_label needs to be provided to take into account the paired structure in cross-validation and permutation testing.")
      }
    }

    # add return values to this list as we go along
    return_values <- list()

    # stores the number of features selected for each fold
    # during cross-validation. we need this later for the
    # random features test
    feats_per_fold <- list()

    # if score is a single function instead of a list of
    # functions, wrap it in a list
    if (class(score) != "list") {
      score <- list(score)
    }
    # ----------------- END INITIAL PROCESSING ----------------- #



    # ----------------- BEGIN CROSS-VALIDATION ----------------- #
    # split data into folds
    if (options$paired) {
      folds <- caret::createFolds(seq(1, nrow(X)/2), options$n_folds)
      fold_names <- names(folds)

      for (fname in fold_names) {
        folds[[fname]] <- which(options$X_label %in% options$X_label[folds[[fname]]])
      }

    } else {
      folds <- caret::createFolds(y, options$n_folds)
      fold_names <- names(folds)
    }

    # vector of cross-validation predictions
    y_pred <- y

    for (fname in fold_names) {
      indices <- folds[[fname]]
      X_train <- X[-indices, , drop = FALSE]
      y_train <- y[-indices]
      X_pred <- X[indices, , drop = FALSE]

      real_features <- select(X_train, y_train)

      # actually, check more for valid indices...
      if (length(real_features) == 0) {
        stop("method$select() did not return any features")
      }

      # store number of features selected in fold for later
      feats_per_fold[[fname]] <- length(real_features)

      model <- train(as.matrix(X_train[, real_features, drop = FALSE]), y_train)
      try_out <- try(y_pred[indices] <- predict(model,
                                                as.matrix(X_pred[, real_features, drop = FALSE])),
                     silent = T)
      if (is(try_out, "try-error")) {
        # No model could be build, provide a model with two latent variables
        # (even if this is not significant)
        opts_model <- list(n_LV = 2)
        model <- train(as.matrix(X_train[, real_features, drop = FALSE]), y_train, opts_model)
        y_pred[indices] <- predict(model,
                                   as.matrix(X_pred[, real_features, drop = FALSE]))
      }
    }
    return_values$cv_y <- y_pred

    # apply the list of functions in score to y_pred, y
    f_star <- function(f) {f(y, y_pred)}
    return_values$cv_score <- lapply(score, f_star)
    print("end CV")
    # ----------------- END CROSS-VALIDATION ----------------- #



    # ----------------- BEGIN RANDOM FEATURES ----------------- #
    n_trials_r <- options$rf_trials
    if (n_trials_r > 0) {
      n_scores <- length(score)
      rf_scores <- list(vector(mode = "numeric", length = n_trials_r))
      rf_scores <- rep(rf_scores, n_scores)

      for (trial_r in 1:n_trials_r) {

        for (fname in fold_names) {
          indices <- folds[[fname]]
          X_train <- X[-indices, , drop = FALSE]
          y_train <- y[-indices]
          X_pred <- X[indices, , drop = FALSE]

          # careful with sample() pathology here...
          # select random features that are NOT ones already chosen
          total_features <- colnames(X)
          nonoverlap_features <- total_features[-which(total_features %in% real_features)]
          ##FEATS PER FOLD IS USING TRIAL_R FNAME AND NOT OVERALL FNAME (wait this is okay)
          random_features <- sample(nonoverlap_features, feats_per_fold[[fname]])
          opts_rf <- list(n_LV = 2)
          model <- train(as.matrix(X_train[, random_features, drop = FALSE]), y_train, opts_rf)

          y_pred[indices] <- predict(model, as.matrix(X_pred[, random_features, drop = FALSE]))
        }

        # compute list of scores
        score_list <- lapply(score, f_star)

        # assign them to vectors in the list
        rf_scores <- lv_assign(rf_scores, score_list, trial_r)
      }

      return_values$rf_scores <- rf_scores
    }
    print("end random features")
    # ----------------- END RANDOM FEATURES ----------------- #



    # ----------------- BEGIN PERMUTATION TESTING ----------------- #
    n_trials_p <- options$pt_trials

    if (n_trials_p > 0) {
      n_scores <- length(score)
      pt_scores <- list(vector(mode = "numeric", length = n_trials_p))
      pt_scores <- rep(pt_scores, n_scores)

      y_perm <- y

      for (trial_p in 1:n_trials_p) {
        if (options$paired) { # create permuted y, flip pairs with a 50% probability
          for (fname in fold_names) {
            indices <- folds[[fname]]
            tmp_rn <- runif(length(indices)/2)
            flip_pairs <- which(tmp_rn > 0.5)
            flip_pairs_labels <- unique(options$X_label[indices])[flip_pairs]

            y_perm[indices] <- y[indices]

            for (ind_pair in flip_pairs_labels) {
              y_perm[which(options$X_label == ind_pair)[2]] <- y[which(options$X_label == ind_pair)[1]]
              y_perm[which(options$X_label == ind_pair)[1]] <- y[which(options$X_label == ind_pair)[2]]
            }
          }
        } else {# create permuted y, but only permute inside each fold
          for (fname in fold_names) {
            indices <- folds[[fname]]
            perm <- sample(1:length(indices))
            y_perm[indices] <- y[indices[perm]]
          }
        }

        for (fname in fold_names) {
          indices <- folds[[fname]]
          X_train <- X[-indices, , drop = FALSE]
          y_train <- y_perm[-indices]
          X_pred <- X[indices, , drop = FALSE]

          features <- select(X_train, y_train)
          model <- train(as.matrix(X_train[, features, drop = FALSE]), y_train)

          try_out <- try(y_pred[indices] <- predict(model,
                                                    as.matrix(X_pred[, features, drop = FALSE])),
                         silent = T)
          if (is(try_out, "try-error")) {
            # No model could be build, provide a model with two latent variables
            # (even if this is not significant)
            opts_model <- list(n_LV = 2)
            model <- train(as.matrix(X_train[, features, drop = FALSE]), y_train, opts_model)
            y_pred[indices] <- predict(model,
                                       as.matrix(X_pred[, features, drop = FALSE]))
          }
        }

        if (options$compare_pred == "y") {
          # compute list of scores
          score_list <- lapply(score, f_star)
        } else if (options$compare_pred == "y_perm") {
          f_star_perm <- function(f) {f(y_perm, y_pred)}
          score_list <- lapply(score, f_star_perm)
        }

        pt_scores <- lv_assign(pt_scores, score_list, trial_p)
      }

      return_values$pt_scores <- pt_scores
    }
    # ----------------- END PERMUTATION TESTING ----------------- #
    print("end permutation testing")


    # ----------------- BEGIN FINAL PROCESSING ----------------- #
    # unpack the list if its length is 1 (score was just a function)
    if (length(score) == 1) {
      return_values$cv_score <- return_values$cv_score[[1]]

      if ("rf_scores" %in% names(return_values)) {
        return_values$rf_scores <- return_values$rf_scores[[1]]
      }

      if ("pt_scores" %in% names(return_values)) {
        return_values$pt_scores <- return_values$pt_scores[[1]]
      }
    }
    # ----------------- END FINAL PROCESSING ----------------- #

    # remove the actual prediction from the validation
    return_values["cv_y"] <- NULL
    vals_list[[trial]] <- return_values
  }
  return(vals_list)
}

# for each index in vec_list, it sets
# vec_list[[ind]][v_index] = val_list[[ind]]
lv_assign <- function(vec_list, val_list, v_index) {
  list_indices <- c(1:length(vec_list))
  for (ind in list_indices) {
    vec <- vec_list[[ind]]
    vec[v_index] <- val_list[[ind]]
    vec_list[[ind]] <- vec
  }
  return(vec_list)
}

#' Visualization of validation results
#'
#' @param vals
#' @param options
#'
#' @return plot handle
#' @export
visualize_validate_specialized <- function(vals, full_val, strain_val, options = list()) {

  # ----------------- OPTIONS ----------------- #
  if (!("y_label" %in% names(options))) {
    options$y_label <- "score"
  }
  # ----------------- END OPTIONS ----------------- #


  tmp_vals <-  unlist(vals)
  tmp_vals <- tmp_vals[which(grepl("score", names(tmp_vals)))]
  df_val <- data.frame(score = tmp_vals,
                       model = gsub("_.*", "", names(tmp_vals)))
  df_val$model <- factor(df_val$model, levels = unique(df_val$model))
  additional <- data.frame(score = c(full_val, unlist(strain_val)),
                           model = c("full model",
                                     rep("leave-one-strain-out",5)))
  df_val <- rbind(df_val, additional)

  # assign x-labels
  x_labels <- levels(df_val$model)
  x_labels[which(x_labels == "cv")] <- "model"
  x_labels[which(x_labels == "rf")] <- "random \nfeatures"
  x_labels[which(x_labels == "pt")] <- "permuted \nlabels"


  # Calculate p-values and generate label
  if (is.null(names(vals))) {
    n_repl <- length(vals)
  } else {
    n_repl <- 1
  }

  if (n_repl > 1) {
    # random features
    if ("rf" %in% levels(df_val$model)) {
      pval_rf <- rep(NA, length = length(vals))

      for (ind in 1:length(vals)) {
        pval_rf[ind] <- length(which(vals[[ind]]$rf_scores > vals[[ind]]$cv_score))/length(vals[[ind]]$rf_scores)
      }

      if (median(pval_rf) == 0) {
        label_rf <- paste0("p<", 1/length(vals[[ind]]$rf_scores))
      } else {
        label_rf <- paste0("p=", mean(pval_rf))
      }
    }

    # permuted labels
    if ("pt" %in% levels(df_val$model)) {
      pval_pt <- rep(NA, length = length(vals))

      for (ind in 1:length(vals)) {
        pval_pt[ind] <- length(which(vals[[ind]]$pt_scores > vals[[ind]]$cv_score))/length(vals[[ind]]$pt_scores)
      }

      if (median(pval_pt) == 0) {
        label_pt <- paste0("p<", 1/length(vals[[ind]]$pt_scores))
      } else {
        label_pt <-  paste0("p=", mean(pval_pt))
      }
    }
  } else {
    # random features
    if ("rf" %in% levels(df_val$model)) {
      pval_rf <- length(which(vals$rf_scores > vals$cv_score))/length(vals$rf_scores)
      if (pval_rf == 0) {
        label_rf <- paste0("p<", 1/length(vals$rf_scores))
      } else {
        label_rf <- paste0("p=", pval_rf)
      }
    }

    # permuted labels
    if ("pt" %in% levels(df_val$model)) {
      pval_pt <- length(which(vals$pt_scores > vals$cv_score))/length(vals$pt_scores)
      if (pval_pt == 0) {
        label_pt <- paste0("p<", 1/length(vals$pt_scores))
      } else {
        label_pt <-  paste0("p=", pval_pt)
      }
    }
  }


  y_pos <- max(df_val$score) + 0.05

  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m - sd(x)
    ymax <- m + sd(x)
    return(c(y = m, ymin = ymin, ymax = ymax))
  }

  plt <- ggplot2::ggplot(df_val, ggplot2::aes(x = model, y = score), fill = "gray") +
    ggplot2::geom_violin(fill = "gray", color = "gray") +
    ggplot2::stat_summary(fun.data = data_summary, geom = "pointrange", size = 0.6, fatten = .8, color = "black") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(color = "black", size = 6),
                   axis.title = ggplot2::element_text(color = "black", size = 8),
                   axis.text.x =  ggplot2::element_text(size = 8, angle = 0, hjust = 0.5)) +
    ggplot2::ylab(options$y_label) +
    ggplot2::scale_x_discrete("", labels = x_labels)

  if ("rf" %in% levels(df_val$model)) {
    plt <- plt +  ggpubr::geom_bracket(xmin = 1, xmax = which(levels(df_val$model) == "rf"),
                                       inherit.aes = FALSE, label.size = 2.5,
                                       y.position = y_pos, label = label_rf)
  }

  if ("pt" %in% levels(df_val$model)) {
    plt <- plt +  ggpubr::geom_bracket(xmin = 1, xmax = which(levels(df_val$model) == "pt"),
                                       inherit.aes = FALSE, label.size = 2.5,
                                       y.position = y_pos + 0.12, label = label_pt)
  }

  if (grepl("ccuracy", options$y_label)) {
    plt <- plt + ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1),
                                             labels = c("0", "0.5", "1"),
                                             limits = c(0, max(1, max(df_val$score) + 0.22)))
  }
  return(plt)

}
