## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE,
                      fig.align = "center")

## ----timeline, fig.width=7, fig.height=3.8, out.width="\\textwidth", fig.cap="glmnet release milestones. Top row runs left-to-right; the path folds down at the right and the bottom row runs right-to-left, ending at v5.0. Evenly spaced for legibility rather than to reflect calendar time."----
par(mar = c(0.4, 0.4, 0.4, 0.4))
milestones <- data.frame(
  ver = c("Initial release (2008)", "1.2 (2010)", "1.6 (2011)",
          "2.0-1 (2015)", "3.0 (2019)",
          "4.0 (2020)", "4.1 (2021)", "4.1-3 (2021)", "5.0 (2026)"),
  lbl = c("coordinate descent;\nFortran kernel\n(Friedman, Hastie,\nTibshirani)",
          "coxnet\n(Simon)",
          "strong rules screen\n(Tibshirani et al.)",
          "CV rewrite\n(Hastie)",
          "relax; assess;\nbigGlm\n(Hastie,\nNarasimhan)",
          "GLM extension\n(Tay, Hastie,\nNarasimhan)",
          "(start, stop] Cox;\nstrata; sparse Cox\n(Tay, Narasimhan,\nHastie)",
          "Fortran to C++\n(Yang)",
          "coxdev submodules;\nCox unified\n(Taylor,\nNarasimhan, Hastie)")
)
n <- nrow(milestones); ntop <- 5; nbot <- n - ntop
top_y <- 1.4; bot_y <- -1.4
col_x <- seq_len(ntop)
px <- numeric(n); py <- numeric(n)
px[1:ntop]     <- col_x;                   py[1:ntop]     <- top_y
px[(ntop+1):n] <- rev(col_x[(ntop-nbot+1):ntop]); py[(ntop+1):n] <- bot_y

plot.new(); plot.window(xlim = c(0.4, ntop + 0.45), ylim = c(-2.3, 2.3))

seg <- function(x0, y0, x1, y1)
  segments(x0, y0, x1, y1, lwd = 1.4, col = "grey35")
arr <- function(x0, y0, x1, y1)
  arrows(x0, y0, x1, y1, length = 0.09, angle = 25, lwd = 1.4,
         col = "grey35")
pad <- 0.18
# top row: left -> right
for (i in 1:(ntop-1)) arr(px[i] + pad, py[i], px[i+1] - pad, py[i+1])
# fold down at right: tiny side step to clear labels on the shared column
fold_x <- ntop + 0.3
seg(px[ntop] + pad, py[ntop], fold_x, py[ntop])        # right from last top node
seg(fold_x, py[ntop], fold_x, py[ntop+1])              # down along margin
arr(fold_x, py[ntop+1], px[ntop+1] + pad, py[ntop+1])  # left into first bottom node
# bottom row: right -> left
for (i in (ntop+1):(n-1)) arr(px[i] - pad, py[i], px[i+1] + pad, py[i+1])

points(px, py, pch = 21, bg = "grey90", col = "grey20",
       cex = 2.0, lwd = 1.3)
text(px, py, labels = seq_len(n), cex = 0.72, col = "grey15",
     font = 2)

# version labels: above top row, below bottom row
text(px[1:ntop], py[1:ntop] + 0.28, milestones$ver[1:ntop],
     cex = 0.72, col = "grey25", pos = 3, offset = 0)
text(px[(ntop+1):n], py[(ntop+1):n] - 0.28, milestones$ver[(ntop+1):n],
     cex = 0.72, col = "grey25", pos = 1, offset = 0)

# description labels: below top row, above bottom row (interior)
text(px[1:ntop], py[1:ntop] - 0.28, milestones$lbl[1:ntop],
     cex = 0.66, pos = 1, offset = 0)
text(px[(ntop+1):n], py[(ntop+1):n] + 0.28, milestones$lbl[(ntop+1):n],
     cex = 0.66, pos = 3, offset = 0)

