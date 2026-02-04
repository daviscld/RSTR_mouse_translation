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

