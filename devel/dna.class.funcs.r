create.seq.df <- function(g, label = NULL) {
  library(strataG)
  stopifnot.gtypes(g, "haploid")
  stopifnot.aligned(g$sequences)
  options(stringsAsFactors = F)

  # extract strata data frame
  g <- decode(g)
  strata.df <- as.data.frame(g$genotypes)
  strata.df$ids <- rownames(strata.df)

  # remove rows without stratifications
  strata.df <- na.omit(strata.df)

  # sequences to use are those in the updated strata data.frame
  dna.seqs <- g$sequences[which(names(g$sequences) %in% strata.df$haplotype)]

  # extract variable sites for these sequences and create sequence matrix
  var.sites <- variable.sites(dna.seqs)
  var.seq.mat <- tolower(do.call(rbind, var.sites$sites))
  sites <- paste("site", colnames(var.sites$site.freqs), sep = ".")

  # only use sites with valid bases
  to.keep <- apply(var.seq.mat, 2, function(x) all(x %in% c("a", "c", "g", "t", "-", ".")))
  if(sum(to.keep) == 0) return(NULL)
  var.seq.mat <- cbind(var.seq.mat[, to.keep])
  colnames(var.seq.mat) <- sites[to.keep]

  # create factors of variable site columns
  var.seq.df <- do.call(data.frame, lapply(colnames(var.seq.mat), function(x) factor(var.seq.mat[, x])))
  colnames(var.seq.df) <- colnames(var.seq.mat)
  var.seq.df <- cbind(haplotype = names(var.sites$sites), var.seq.df)
  seq.df <- merge(strata.df, var.seq.df, all = T, by = "haplotype")

  # remove any rows with missing data and remove id and haplotype columns
  seq.df <- na.omit(seq.df)
  seq.df$strata <- factor(seq.df$strata)
  rownames(seq.df) <- seq.df$ids
  seq.df$id <- seq.df$haplotype <- NULL
  sites <- seq.df[, -1, drop = FALSE]

  # remove sites where substitutions are only represented by one individual
  to.keep <- apply(sites, 2, function(x) sum(table(x) > 1) > 1)
  if(sum(to.keep) == 0) return(NULL)

  strata <- if(is.null(label)) seq.df$strata else paste(label, seq.df$strata)
  cbind(strata = factor(strata), sites[, to.keep, drop = FALSE])
}

classify.by.hap.freq <- function(x) {
  x.haps <- make.haps(x[, -1, drop = F])
  hap.freq <- prop.table(table(x.haps, x$strata), 1)
  pred <- sapply(x.haps, function(i) colnames(hap.freq)[which.max(hap.freq[as.character(i), ])])
  strata.levels <- unique(x$strata)
  class.tbl <- table(factor(x$strata, levels = strata.levels), factor(pred, levels = strata.levels))
  overall.diag <- sum(diag(class.tbl)) / sum(class.tbl)
  class.tbl <- cbind(class.tbl, diagnosability = diag(class.tbl) / rowSums(class.tbl))
  list(pred = pred, class.tbl = class.tbl, overall.diag = overall.diag)
}


calc.minimum.error <- function(x) {
  x.tbl <- table(make.haps(x[, -1, drop = F]), x[, 1])
  min.freq.shared.haps <- rowSums(apply(x.tbl, 1, function(y) {
    y[which.max(y)] <- 0
    y
  }))
  sum(min.freq.shared.haps) / sum(x.tbl)
}

make.haps <- function(x) {
  dna.seq <- apply(x, 1, paste, collapse = "")
  most.freq <- sort(table(dna.seq), decreasing = T)
  as.numeric(factor(dna.seq, levels = names(most.freq)))
}

collapse.to.haps <- function(x) {
  x.haps <- make.haps(x)
  x[!duplicated(x.haps), , drop = F]
}

hap.info <- function(x) {
  x.haps <- make.haps(x[, -1])
  c(num.haps = max(x.haps), div = calc.diversity(x.haps))
}

rf.test <- function(x, ...) {
  stopifnot(require(randomForest))
  stopifnot(require(reshape2))
  stopifnot(require(ggplot2))
  rf <- randomForest(strata ~ ., data = x, ...)
  rf.votes <- cbind(strata = x$strata, data.frame(rf$votes, check.names = FALSE))
  rf.votes <- melt(rf.votes, id.vars = "strata", variable.name = "predicted", value.name = "votes")
  rf.votes <- rf.votes[order(rf.votes$strata, rf.votes$predicted, rf.votes$votes), ]
  p <- if(nlevels(x$strata) > 2) {
    ggplot(rf.votes, aes(predicted, votes)) + geom_violin() + geom_jitter() + facet_grid(~ strata)
  } else {
    rf.votes <- subset(rf.votes, as.character(strata) == as.character(predicted))
    ggplot(rf.votes, aes(strata, votes)) + geom_violin() + geom_jitter()
  }
  print(p)
  varImpPlot(rf)
  print(rf)
  invisible(rf)
}

min.votes <- function(rf, mv.vec) {
  vote.dfs <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    df <- data.frame(votes = votes, is.correct = rf$predicted[y.i] == rf$y[y.i])
    df[order(df$votes), ]
  }, simplify = F)
  sapply(1 - mv.vec, function(i) {
    min(sapply(vote.dfs, function(df) {
      cum.pct <- 1:nrow(df) / nrow(df)
      df$votes[min(which(cum.pct >= i))]
    }))
  })
}

pct.diag <- function(rf, pd.vec) {
  vote.dfs <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    df <- data.frame(votes = votes, is.correct = rf$predicted[y.i] == rf$y[y.i])
    df[order(df$votes), ]
  }, simplify = F)
  strata.pd <- sapply(vote.dfs, function(df) sum(df$is.correct) / nrow(df))
  is.min <- which.min(strata.pd)
  min.df <- vote.dfs[[is.min]]
  result <- sapply(pd.vec, function(i) sum(min.df$votes > i) / nrow(min.df))
  result <- c(strata.pd[is.min], result)
  names(result) <- paste("pd", c("0.5", pd.vec), sep = ".")
  result
}

plot.vote.dist <- function(rf, col) {
  obs.votes <- sapply(colnames(rf$votes), function(y) {
    y.i <- which(rf$y == y)
    votes <- rf$votes[y.i, y]
    df <- data.frame(votes = votes, is.correct = rf$predicted[y.i] == rf$y[y.i])
    df[order(df$votes), ]
  }, simplify = F)
  vote.vals <- unlist(lapply(obs.votes, function(y) y$votes))
  op <- par(mfrow = c(length(obs.votes), 1), mar = c(3, 3, 2, 1) + 0.1, oma = c(2, 2, 0, 0))
  for(x in names(obs.votes)) {
    vote.thresh <- 1 / length(obs.votes)
    votes <- obs.votes[[x]]$votes
    is.correct <- obs.votes[[x]]$is.correct
    cum.pct <- 1:length(votes) / length(votes)
    plot(cum.pct, votes, type = "l", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
    abline(v = 0.05, lwd = 2, lty = "dashed")
    abline(v = 0.25, lwd = 2, lty = "dashed")
    abline(h = 0.95, lwd = 2, lty = "dashed")
    points(cum.pct[is.correct], votes[is.correct], pch = 21, col = "red", bg = "red")
    points(cum.pct[!is.correct], votes[!is.correct], pch = 21, col = "black", bg = "black")
    mtext(x, side = 3, line = 0.5, adj = 0, col = col)
  }
  mtext("Cumulative % of samples", side = 1, line = 0.5, outer = T)
  mtext("% Votes", side = 2, line = 0.5, outer = T)
  par(op)
}

unique.hap.by.site <- function(x) {
  x.haps <- make.haps(x)
  result <- unlist(apply(x[, -1, drop = F], 2, function(b) {
    site.tbl <- table(x.haps, b)
    shared.site <- apply(site.tbl, 2, function(y) sum(y > 0))
    j <- which(shared.site == 1)
    if(length(j) > 0) {
      rownames(site.tbl)[which(site.tbl[, j] > 0)]
    } else NA
  }))
  result[!is.na(result)]
}


calc.conf.ci <- function(rf) {
  conf <- rf$confusion
  t(apply(conf, 1, function(y) {
    diagnosable <- as.vector(1 - y[length(y)])
    y <- y[-length(y)]
    ci <- binom.test(sum(y) * diagnosable, sum(y), diagnosable)$conf.int
    names(ci) <- c("LCI", "UCI")
    c(diagnosable = diagnosable, ci)
  }))
}

genetic.rf <- function(seq.df, ...) {
  if(is.null(seq.df)) return(NULL)
  stopifnot(require(randomForest))
  rf <- randomForest(strata ~ ., data = seq.df, ...)
  
  overall.accuracy <- 1 - as.vector(rf$err.rate[nrow(rf$err.rate), "OOB"])
  conf.ci <- calc.conf.ci(rf)
  min.diag <- which.min(conf.ci[, 1])
  diag.strata <- rownames(conf.ci)[min.diag]
  prior.diag <- as.vector(prop.table(table(rf$y))[diag.strata])
  hap.class.tbl <- classify.by.hap.freq(seq.df)$class.tbl
  shared.hap.diag <- hap.class.tbl[diag.strata, "diagnosability"]
  haps <- make.haps(seq.df[, -1, drop = FALSE])
  num.haps <- length(unique(haps))
  vs <- ncol(seq.df) - 1
  unique.hap.vs.vec <- unique.hap.by.site(seq.df)
  num.rf.vs.unique.hap <- length(unique.hap.vs.vec)
  num.unique.hap.rf.vs <- length(table(unique.hap.vs.vec))

  smry <- data.frame(
    overall.accuracy = overall.accuracy,
    diag.strata = diag.strata,
    diagnosability = conf.ci[min.diag, 1],
    diag.lci = conf.ci[min.diag, 2],
    diag.uci = conf.ci[min.diag, 3],
    prior.diag = prior.diag,
    shared.hap.diag = shared.hap.diag,
    pd95 = pct.diag(rf, 0.95)[2],
    num.vs = vs,
    num.haps = num.haps,
    hap.div = diversity(haps),
    num.vs.unique.hap = num.rf.vs.unique.hap,
    pct.vs.haps.shared = (vs - num.rf.vs.unique.hap) / vs,
    num.unique.haps.vs = num.unique.hap.rf.vs,
    pct.haps.vs.shared = (num.haps - num.unique.hap.rf.vs) / num.haps
  )
  rownames(smry) <- NULL
  
  list(smry = smry, rf = rf)
}

gtype.rf <- function(g, pairwise = FALSE, ...) {
  rf.func <- function(g) {
    seq.df <- create.seq.df(g)
    strata.freq <- table(seq.df$strata)
    min.n <- min(strata.freq)
    genetic.rf(seq.df, replace = TRUE, sampsize = rep(min.n, length(strata.freq)), ...)
  }
  
  if(pairwise) {
    sp <- strata.pairs(g)
    result <- lapply(1:nrow(sp), function(i) {
      pair.g <- subset.gtypes(g, strata = unlist(sp[i, ]))
      gtype.rf(pair.g, ...)
    })
    list(smry = cbind(sp, do.call(rbind, lapply(result, function(x) x$smry))),
         rf = lapply(result, function(x) x$rf)
    )
  } else rf.func(g)
}

vote.plot <- function(df) {
  df <- cbind(df[, 1], df[, sort(colnames(df)[-1], dec = T)])
  order.list = c(list(df[, 1]), lapply(2:ncol(df), function(i) df[, i]))
  df.order <- do.call(order, order.list)
  df <- df[rev(df.order), ]
  col <- rainbow(length(unique(df[, 1])))
  names(col) <- unique(df[, 1])
  bp <- barplot(t(as.matrix(df[, -1])), horiz = T, col = col, axes = F, axisnames = F, space = 0, border = NA, angle = 0)
  tx <- sort(tapply(bp, df[, 1], min) + 0.1)
  tx <- c(tx - tx[1], max(bp) + tx[1])
  axis(2, at = tx, labels = FALSE, lwd = 2)
  abline(h = tx)
  lbl.x <- sapply(1:(length(tx) - 1), function(i) tx[i] + (tx[i + 1] - tx[i]) / 2)
  names(lbl.x) <- names(tx)[1:(length(tx) - 1)]
  mtext(names(lbl.x), side = 2, at = lbl.x, line = 1, las = 2, col = col[names(lbl.x)], font = 2)
  axis(1, las = 1, lwd = 2)
  mtext("%Votes", side = 1, line = 2, cex = 1)
}

rf.species.id <- function(g, ref.strata, unk.strata, ...) {
  if(length(intersect(ref.strata, unk.strata)) > 0) stop("'ref.strata' and 'unk.strata' cannot overlap.")
  if(length(unique(ref.strata)) < 2) stop("'ref.strata' must have 2 or more strata.")
  sub.g <- subset(g, strata = c(ref.strata, unk.strata))
  seq.df <- create.seq.df(sub.g)
  ref.df <- subset(seq.df, strata %in% ref.strata)
  ref.df$strata <- droplevels(ref.df$strata)
  unk.df <- subset(seq.df, strata %in% unk.strata)
  unk.df$strata <- droplevels(unk.df$strata)
  ref.rf <- genetic.rf(ref.df, ...)
  unk.pred <- predict(ref.rf$rf, unk.df, type = "response")
  unk.prob <- predict(ref.rf$rf, unk.df, type = "prob")
  list(pred = unk.pred, prob = unk.prob, rf = ref.rf)
}
