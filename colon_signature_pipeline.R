###############################################################################
# Colon/rectal cancer signature pipeline (sigident-style; Vey et al. 2019)
#
# Data
#   GSE39582 : colon tumors + 19 non-tumoral mucosa, Affy GPL570, has survival
#   GSE17536 : colorectal tumors only (177), Affy GPL570, has survival
#   GSE87211 : RECTAL tumors (203) + 160 matched mucosa, Agilent GPL13497
#
# Workflow
#   1. Download + phenotype parsing (CHECK the printed tables before going on)
#   2. Probe -> gene symbol, log2 check
#   3. Merge the two GPL570 series, PCA, ComBat batch correction
#   4. Differential expression (limma, BH FDR) + heatmap
#   5. GO / KEGG enrichment (clusterProfiler)
#   6. Diagnostic signature: LASSO + elastic net, 80/20 split, external test
#      in GSE87211 (cross-platform, gene level)
#   7. Prognostic signature: univariate Cox -> penalized Cox risk score trained
#      in GSE39582, validated in GSE17536
#
# TCGA is deliberately left out; add it later as another validation cohort.
###############################################################################

## ---- 0. Setup --------------------------------------------------------------
cran_pkgs <- c("glmnet", "survival", "survminer", "pROC", "caret",
               "pheatmap", "ggplot2")
bioc_pkgs <- c("GEOquery", "Biobase", "limma", "sva", "clusterProfiler",
               "org.Hs.eg.db", "hgu133plus2.db", "AnnotationDbi")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)

# IMPORTANT: if packages were installed or updated in this R session, RESTART R
# (RStudio: Session > Restart R, or Ctrl+Shift+F10) and run the script again from
# a fresh session. Otherwise you get errors like "namespace 'ggplot2' cannot be
# unloaded" or "rlang ... is already loaded, but >= x is required".
# Best practice: run it with source("colon_signature_pipeline.R") so it stops at
# the first error instead of cascading.
# Project folder: the script moves here itself, so no manual setwd() is needed.
# Use forward slashes. Edit this path if you move the project.
PROJECT_DIR <- "C:/Users/alexp/Documents/Research Project/download colon"
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR not found: ", PROJECT_DIR)
setwd(PROJECT_DIR)
message("Working directory: ", getwd())

options(timeout = 1800)   # GEO series matrices are large; default is 60 s

library(ggplot2)          # load first to avoid version clashes with survminer/caret
suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(limma); library(sva)
  library(glmnet); library(survival); library(survminer); library(pROC)
  library(caret); library(pheatmap)
  library(clusterProfiler); library(org.Hs.eg.db)
})

set.seed(2026)
dir.create("results", showWarnings = FALSE)
dir.create("geo_cache", showWarnings = FALSE)

Q_CUT   <- 0.05   # FDR cutoff for DEGs
LFC_CUT <- 1      # |log2FC| cutoff (the lung paper used 2; 1 is more usual)

# Regexes used to find survival columns. Adjust after inspecting the pData
# column names printed in section 1.
SURV_PATTERNS <- list(
  GSE39582 = list(time = "os\\.delay",             event = "os\\.event"),
  GSE17536 = list(time = "overall survival follow", event = "overall_event")
)

## ---- Helpers ---------------------------------------------------------------
# Files are cached in ./geo_cache. If a download keeps failing, download
# <GSE>_series_matrix.txt.gz from the GEO page in your browser, put it in
# ./geo_cache, and re-run: getGEO will use the local copy.
# getGPL is only needed for the Agilent series (gene symbols come from the
# platform table); the Affymetrix ones are annotated with hgu133plus2.db.
get_eset <- function(gse, gpl = FALSE) {
  f <- file.path("geo_cache", paste0(gse, "_series_matrix.txt.gz"))
  if (file.exists(f)) {
    message("Reading local file ", f, " (no network needed)")
    g <- getGEO(filename = f, GSEMatrix = TRUE, getGPL = gpl, destdir = "geo_cache")
  } else {
    g <- getGEO(gse, GSEMatrix = TRUE, getGPL = gpl, destdir = "geo_cache")
  }
  if (is.list(g)) g[[which.max(sapply(g, ncol))]] else g
}

# Tumor / normal labels from the text of the characteristics columns.
assign_group <- function(es) {
  pd   <- pData(es)
  cols <- unique(c(intersect("title", colnames(pd)),
                   grep(":ch1$", colnames(pd), value = TRUE)))
  txt  <- apply(pd[, cols, drop = FALSE], 1, paste, collapse = " | ")
  is_normal <- grepl("non[- ]?tumou?r|mucosa|normal|control|adjacent",
                     txt, ignore.case = TRUE)
  factor(ifelse(is_normal, "normal", "tumor"), levels = c("normal", "tumor"))
}

# Make sure values are on a log2 scale.
ensure_log2 <- function(m) {
  if (quantile(m, 0.99, na.rm = TRUE) > 100) {
    message("  values look un-logged; applying log2(x + 1)")
    m <- log2(pmax(m, 0) + 1)
  }
  m
}

# Collapse probes to one row per gene symbol (probe with highest mean).
collapse_to_gene <- function(mat, genes) {
  keep <- !is.na(genes) & genes != "" & !grepl("///", genes, fixed = TRUE)
  mat <- mat[keep, , drop = FALSE]; genes <- genes[keep]
  ord <- order(rowMeans(mat), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; genes <- genes[ord]
  dup <- duplicated(genes)
  mat <- mat[!dup, , drop = FALSE]
  rownames(mat) <- genes[!dup]
  mat
}

gpl570_to_gene <- function(es) {
  m   <- ensure_log2(exprs(es))
  sym <- AnnotationDbi::mapIds(hgu133plus2.db::hgu133plus2.db,
                               keys = rownames(m), column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
  collapse_to_gene(m, unname(sym))
}

agilent_to_gene <- function(es) {
  m  <- ensure_log2(exprs(es))
  fd <- fData(es)
  sc <- grep("^GENE_SYMBOL$|gene.?symbol", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sc)) stop("No gene symbol column in fData; inspect colnames(fData(es)).")
  collapse_to_gene(m, as.character(fd[[sc]]))
}

find_col <- function(pd, pattern) {
  hit <- grep(pattern, colnames(pd), ignore.case = TRUE, value = TRUE)
  if (!length(hit)) stop("No pData column matches '", pattern, "'. Inspect colnames(pData(es)).")
  hit[1]
}

get_surv <- function(es, patterns) {
  pd <- pData(es)
  tc <- find_col(pd, patterns$time); ec <- find_col(pd, patterns$event)
  message("  survival columns: time = '", tc, "', event = '", ec, "'")
  raw_t <- as.character(pd[[tc]]); raw_e <- as.character(pd[[ec]])
  message("  example time values : ", paste(head(unique(raw_t), 5), collapse = " | "))
  message("  event values        : ",
          paste(names(table(raw_e, useNA = "ifany")), table(raw_e, useNA = "ifany"),
                sep = "=", collapse = ", "))
  # time: pull the first number out of the string
  tm <- suppressWarnings(as.numeric(sub("^[^0-9.-]*(-?[0-9.]+).*$", "\\1", raw_t)))
  # event: accept 0/1 or text such as dead/alive
  ev_txt <- tolower(trimws(raw_e))
  ev <- rep(NA_real_, length(ev_txt))
  ev[grepl("^1|dead|deceas|^yes|^true|death", ev_txt)] <- 1
  # negations override: "no death", "not dead", "alive", "censored", 0 ...
  ev[grepl("^0|alive|^no\\b|^not\\b|none|without|^false|censor", ev_txt)] <- 0
  if (all(is.na(tm)) || all(is.na(ev)))
    stop("Could not parse survival values for ", tc, " / ", ec,
         ". Check the 'example values' printed above and edit get_surv().")
  if (length(unique(na.omit(ev))) < 2)
    stop("All patients have the same event status after parsing (", paste(unique(ev), collapse = ","),
         "). The event labels printed above are not being interpreted correctly; edit get_surv().")
  message("  parsed events: ", sum(ev == 1, na.rm = TRUE), " events / ", sum(!is.na(ev)), " patients")
  data.frame(sample = colnames(es), time = tm, event = ev, stringsAsFactors = FALSE)
}

zfit   <- function(X) list(mu = colMeans(X), sd = apply(X, 2, sd))
zapply <- function(X, z) sweep(sweep(X, 2, z$mu), 2, pmax(z$sd, 1e-6), "/")

run_limma <- function(expr, group, covar = NULL) {
  if (is.null(covar)) {
    design <- model.matrix(~ 0 + group)
    colnames(design) <- c("normal", "tumor")
  } else {
    design <- model.matrix(~ 0 + group + covar)
    colnames(design)[1:2] <- c("normal", "tumor")
  }
  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, makeContrasts(tumor - normal, levels = design))
  fit <- eBayes(fit)
  tt  <- topTable(fit, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt
}

## ---- 1. Download and inspect ----------------------------------------------
es39582 <- get_eset("GSE39582")
es17536 <- get_eset("GSE17536")
es87211 <- get_eset("GSE87211", gpl = TRUE)   # Agilent: need platform table for symbols

for (nm in c("es39582", "es17536", "es87211")) {
  es <- get(nm)
  cat("\n====", nm, "| platform:", annotation(es), "| samples:", ncol(es), "====\n")
  print(colnames(pData(es)))
}

# Platform sanity checks (merge in section 3 requires both be GPL570)
stopifnot(annotation(es39582) == "GPL570", annotation(es17536) == "GPL570")
if (annotation(es87211) != "GPL13497")
  warning("GSE87211 platform is ", annotation(es87211), ", expected GPL13497")

grp39582 <- assign_group(es39582)
grp17536 <- assign_group(es17536)
grp87211 <- assign_group(es87211)

cat("\nGSE39582 (expect ~566 tumor / 19 normal):\n");  print(table(grp39582))
cat("GSE17536 (expect 177 tumor / 0 normal):\n");      print(table(grp17536))
cat("GSE87211 (expect 203 tumor / 160 normal):\n");    print(table(grp87211))

# >>> STOP HERE on first run. If the counts above are wrong, look at the
# >>> characteristics columns and fix assign_group() before continuing.
# e.g. lapply(pData(es39582)[, grep(":ch1$", colnames(pData(es39582)))], function(x) head(unique(x)))

## ---- 2. Probe -> gene ------------------------------------------------------
g39582 <- gpl570_to_gene(es39582)
g17536 <- gpl570_to_gene(es17536)
g87211 <- agilent_to_gene(es87211)
cat("\nGenes:", nrow(g39582), "/", nrow(g17536), "/", nrow(g87211), "\n")

## ---- 3. Merge GPL570 series, batch-correct ---------------------------------
common <- intersect(rownames(g39582), rownames(g17536))
expr   <- cbind(g39582[common, ], g17536[common, ])
group  <- factor(c(as.character(grp39582), as.character(grp17536)),
                 levels = c("normal", "tumor"))
batch  <- factor(rep(c("GSE39582", "GSE17536"), c(ncol(g39582), ncol(g17536))))
names(group) <- names(batch) <- colnames(expr)

# drop near-constant genes (ComBat fails on zero-variance genes within a batch)
ok <- apply(expr, 1, function(x) all(tapply(x, batch, var) > 1e-8))
expr <- expr[ok, ]
cat("Merged matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

pca_plot <- function(m, title) {
  p <- prcomp(t(m), scale. = TRUE)
  df <- data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], batch = batch, group = group)
  ggplot(df, aes(PC1, PC2, colour = batch, shape = group)) + geom_point(alpha = .7) +
    ggtitle(title) + theme_bw()
}
ggsave("results/pca_before_combat.png", pca_plot(expr, "Before ComBat"), width = 6, height = 5)

# Optional guided PCA batch test, as in the paper (skipped if gPCA is missing)
if (requireNamespace("gPCA", quietly = TRUE)) {
  gp <- gPCA::gPCA.batchdetect(x = t(expr), batch = as.numeric(batch), center = TRUE)
  cat("gPCA delta:", gp$delta, " p:", gp$p.val, "\n")
}

# ComBat with the tumor/normal label protected. NOTE: GSE17536 has no normals,
# so batch and group are partly confounded: batch effects are estimated mainly
# from tumor samples. Interpret DEGs with that in mind.
expr_bc <- ComBat(dat = expr, batch = batch, mod = model.matrix(~ group))
ggsave("results/pca_after_combat.png", pca_plot(expr_bc, "After ComBat"), width = 6, height = 5)

## ---- 4. Differential expression -------------------------------------------
deg_all <- run_limma(expr_bc, group)
deg_sig <- subset(deg_all, adj.P.Val < Q_CUT & abs(logFC) > LFC_CUT)
cat("\nDEGs (q <", Q_CUT, ", |log2FC| >", LFC_CUT, "):", nrow(deg_sig),
    "(", sum(deg_sig$logFC > 0), "up /", sum(deg_sig$logFC < 0), "down )\n")
write.csv(deg_all, "results/DEG_all_genes_merged_GPL570.csv", row.names = FALSE)
write.csv(deg_sig, "results/DEG_significant_merged_GPL570.csv", row.names = FALSE)

# Sensitivity check: batch as a covariate instead of ComBat
deg_cov <- run_limma(expr, group, covar = batch)
sig_cov <- subset(deg_cov, adj.P.Val < Q_CUT & abs(logFC) > LFC_CUT)
cat("Overlap with batch-covariate analysis:",
    length(intersect(deg_sig$gene, sig_cov$gene)), "of", nrow(deg_sig), "\n")

# Heatmap (top 200 by |logFC| if there are many)
hm_genes <- head(deg_sig$gene[order(-abs(deg_sig$logFC))], 200)
if (length(hm_genes) > 1) {
  ann <- data.frame(group = group, batch = batch, row.names = colnames(expr_bc))
  dir.create("results", showWarnings = FALSE)
  tryCatch(
    pheatmap(expr_bc[hm_genes, ], scale = "row", show_colnames = FALSE,
             show_rownames = FALSE, annotation_col = ann,
             filename = "results/heatmap_DEGs.png", width = 8, height = 8),
    error = function(e) message("Heatmap not saved (", conditionMessage(e),
                                "). Check getwd(), and that the file isn't open elsewhere.")
  )
}

# Independent differential expression in GSE87211 (rectal, Agilent)
deg87 <- run_limma(g87211, grp87211)
shared <- intersect(deg_all$gene, deg87$gene)
lfc_a  <- deg_all$logFC[match(shared, deg_all$gene)]
lfc_b  <- deg87$logFC[match(shared, deg87$gene)]
cat("\nlogFC correlation, merged GPL570 vs GSE87211:",
    round(cor(lfc_a, lfc_b, use = "complete.obs"), 3), "\n")
sig_both <- intersect(deg_sig$gene, subset(deg87, adj.P.Val < Q_CUT & abs(logFC) > LFC_CUT)$gene)
sign_ok <- sign(deg_all$logFC[match(sig_both, deg_all$gene)]) ==
           sign(deg87$logFC[match(sig_both, deg87$gene)])
cat("Genes significant in both:", length(sig_both),
    "| same direction:", sum(sign_ok), "\n")
write.csv(deg87, "results/DEG_all_genes_GSE87211.csv", row.names = FALSE)

## ---- 5. Enrichment ---------------------------------------------------------
universe <- bitr(rownames(expr_bc), "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
deg_ids  <- bitr(deg_sig$gene,      "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID

ego <- enrichGO(deg_ids, OrgDb = org.Hs.eg.db, ont = "BP", universe = universe,
                pAdjustMethod = "BH", qvalueCutoff = Q_CUT, readable = TRUE)
ekg <- tryCatch(enrichKEGG(deg_ids, organism = "hsa", universe = universe,
                           pAdjustMethod = "BH", qvalueCutoff = Q_CUT),
                error = function(e) { message("KEGG failed (needs internet): ", e$message); NULL })
write.csv(as.data.frame(ego), "results/GO_BP_enrichment.csv", row.names = FALSE)
if (!is.null(ekg)) write.csv(as.data.frame(ekg), "results/KEGG_enrichment.csv", row.names = FALSE)
ggsave("results/GO_dotplot.png", dotplot(ego, showCategory = 15), width = 8, height = 7)

## ---- 6. Diagnostic signature (LASSO + elastic net) -------------------------
# Features limited to genes present on BOTH platforms so the model can be
# applied to GSE87211.
feat <- intersect(rownames(expr_bc), rownames(g87211))
X <- t(expr_bc[feat, ]); y <- group

strata <- droplevels(interaction(batch, y))
tr <- as.vector(createDataPartition(strata, p = 0.8, list = FALSE))
Xtr <- X[tr, ];  ytr <- y[tr]
Xte <- X[-tr, ]; yte <- y[-tr]
cat("\nTrain:", table(ytr), " Test:", table(yte), "(normal, tumor)\n")

ztr  <- zfit(Xtr)
Xtrz <- zapply(Xtr, ztr); Xtez <- zapply(Xte, ztr)

# External set: z-score within the cohort itself (platform differs)
Xval  <- t(g87211[feat, ])
Xvalz <- zapply(Xval, zfit(Xval)); yval <- grp87211

# Class weights: only ~19 normals vs ~700 tumors
w <- ifelse(ytr == "normal", sum(ytr == "tumor") / sum(ytr == "normal"), 1)

fit_glm <- function(alpha) {
  foldid <- createFolds(ytr, k = 10, list = FALSE)
  cv.glmnet(Xtrz, ytr, family = "binomial", alpha = alpha, weights = w,
            foldid = foldid, type.measure = "deviance")
}

eval_glm <- function(cvfit, s, Xnew, ynew, thr = 0.5) {
  prob <- as.numeric(predict(cvfit, Xnew, s = s, type = "response"))
  r    <- roc(ynew, prob, levels = c("normal", "tumor"), direction = "<", quiet = TRUE)
  pred <- factor(ifelse(prob > thr, "tumor", "normal"), levels = levels(ynew))
  list(cm = confusionMatrix(pred, ynew, positive = "tumor"),
       auc = as.numeric(auc(r)), roc = r, prob = prob)
}

sig_genes <- function(cvfit, s) {
  b <- as.matrix(coef(cvfit, s = s)); rownames(b)[b[, 1] != 0 & rownames(b) != "(Intercept)"]
}

models <- list(LASSO = fit_glm(1), ElasticNet_a0.9 = fit_glm(0.9))
diag_sigs <- list()
for (nm in names(models)) {
  for (s in c("lambda.min", "lambda.1se")) {
    key <- paste(nm, s, sep = "_")
    diag_sigs[[key]] <- sig_genes(models[[nm]], s)
    te  <- eval_glm(models[[nm]], s, Xtez,  yte)
    ext <- eval_glm(models[[nm]], s, Xvalz, yval)
    # Cutoff tuned on GSE87211 itself (OPTIMISTIC; shown because within-cohort
    # z-scoring shifts the probabilities, so 0.5 is not meaningful there)
    best <- as.numeric(coords(ext$roc, "best", best.method = "youden",
                              ret = "threshold", transpose = TRUE))[1]
    ext_y <- eval_glm(models[[nm]], s, Xvalz, yval, thr = best)
    cat(sprintf("\n[%s] %d genes | internal test AUC %.3f, bal.acc %.3f | GSE87211 AUC %.3f\n",
                key, length(diag_sigs[[key]]), te$auc,
                te$cm$byClass["Balanced Accuracy"], ext$auc))
    cat(sprintf("   GSE87211 bal.acc @0.5 = %.3f | @Youden cutoff %.3f = %.3f (optimistic)\n",
                ext$cm$byClass["Balanced Accuracy"], best,
                ext_y$cm$byClass["Balanced Accuracy"]))
    print(ext$cm$table)
  }
}
cat("\nLASSO vs elastic net (lambda.min) overlap:",
    length(intersect(diag_sigs$LASSO_lambda.min, diag_sigs$ElasticNet_a0.9_lambda.min)), "genes\n")

sig_df <- data.frame(model = rep(names(diag_sigs), lengths(diag_sigs)),
                     gene = unlist(diag_sigs), row.names = NULL)
write.csv(sig_df, "results/diagnostic_signatures.csv", row.names = FALSE)
saveRDS(models, "results/diagnostic_glmnet_models.rds")

## ---- 7. Prognostic signature ----------------------------------------------
# Train in GSE39582 (many events), validate in GSE17536.
s39582 <- get_surv(es39582, SURV_PATTERNS$GSE39582)
s17536 <- get_surv(es17536, SURV_PATTERNS$GSE17536)

prep_cohort <- function(sdf, gmat, grp) {
  ok <- grp == "tumor" & !is.na(sdf$time) & !is.na(sdf$event) & sdf$time > 0
  if (sum(ok) == 0)
    stop("No usable survival samples (tumor with time > 0 and non-missing event). ",
         "Check the parsed time/event values printed by get_surv().")
  sdf <- sdf[ok, ]
  X   <- zapply(t(gmat[, sdf$sample]), zfit(t(gmat[, sdf$sample])))  # z-score within cohort
  list(df = sdf, X = X)
}
trn <- prep_cohort(s39582, expr_bc[, colnames(g39582)], grp39582)
val <- prep_cohort(s17536, expr_bc[, colnames(g17536)], grp17536)
cat("\nSurvival train (GSE39582): n =", nrow(trn$df), "events =", sum(trn$df$event),
    "\nSurvival valid (GSE17536): n =", nrow(val$df), "events =", sum(val$df$event), "\n")

# 7a. Univariate Cox on the DEGs (as in the paper), p < 0.05
cand <- intersect(deg_sig$gene, colnames(trn$X))
uni <- t(vapply(cand, function(g) {
  f <- coxph(Surv(trn$df$time, trn$df$event) ~ trn$X[, g])
  s <- summary(f)$coefficients
  c(HR = s[1, "exp(coef)"], p = s[1, "Pr(>|z|)"])
}, numeric(2)))
uni <- data.frame(gene = cand, uni, padj = p.adjust(uni[, "p"], "BH"))
uni_sig <- subset(uni, p < 0.05)
cat("Univariate Cox hits (p < 0.05):", nrow(uni_sig),
    "| BH q < 0.05:", sum(uni$padj < 0.05), "\n")
write.csv(uni[order(uni$p), ], "results/univariate_cox_GSE39582.csv", row.names = FALSE)

# 7b. Penalized Cox risk score (replaces the paper's tumor-vs-normal classifier)
ysurv <- Surv(trn$df$time, trn$df$event)
cvcox <- cv.glmnet(trn$X[, uni_sig$gene], ysurv, family = "cox", alpha = 0.9, nfolds = 10)
s_cox <- "lambda.1se"
b <- as.matrix(coef(cvcox, s = s_cox))
if (sum(b[, 1] != 0) == 0) {
  message("lambda.1se selected no genes; falling back to lambda.min")
  s_cox <- "lambda.min"
  b <- as.matrix(coef(cvcox, s = s_cox))
}
prog_genes <- rownames(b)[b[, 1] != 0]
if (length(prog_genes) == 0)
  stop("Penalized Cox selected no genes even at lambda.min: the candidate genes carry ",
       "no cross-validated prognostic signal in this cohort.")
cat("Prognostic signature (", s_cox, "):", length(prog_genes), "genes\n")
print(round(b[prog_genes, , drop = FALSE], 3))
write.csv(data.frame(gene = prog_genes, coef = b[prog_genes, 1]),
          "results/prognostic_signature.csv", row.names = FALSE)

risk <- function(X) as.numeric(predict(cvcox, X[, uni_sig$gene], s = s_cox, type = "link"))
trn$df$risk <- risk(trn$X); val$df$risk <- risk(val$X)
cutoff <- median(trn$df$risk)                      # threshold fixed on training data
trn$df$grp <- factor(ifelse(trn$df$risk > cutoff, "high", "low"), levels = c("low", "high"))
val$df$grp <- factor(ifelse(val$df$risk > cutoff, "high", "low"), levels = c("low", "high"))

report <- function(d, label) {
  lr <- survdiff(Surv(time, event) ~ grp, data = d)
  p  <- 1 - pchisq(lr$chisq, 1)
  hr <- summary(coxph(Surv(time, event) ~ grp, data = d))$conf.int["grphigh", ]
  cat(sprintf("%s: log-rank p = %.4g | HR(high vs low) = %.2f [%.2f-%.2f]\n",
              label, p, hr["exp(coef)"], hr["lower .95"], hr["upper .95"]))
  # Base-R KM plot (no survminer/xfun dependency)
  tryCatch({
    png(sprintf("results/KM_%s.png", label), width = 900, height = 700, res = 120)
    on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
    fit <- survfit(Surv(time, event) ~ grp, data = d)
    plot(fit, col = c("steelblue", "firebrick"), lwd = 2, mark.time = TRUE,
         xlab = "Months", ylab = "Survival probability", main = label)
    legend("bottomleft", legend = paste0(c("low", "high"), " risk (n=", fit$n, ")"),
           col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
    legend("topright", legend = sprintf("log-rank p = %.3g\nHR = %.2f", p, hr["exp(coef)"]),
           bty = "n")
  }, error = function(e) message("KM plot not saved for ", label, ": ", conditionMessage(e)))
}
report(trn$df, "GSE39582_train"); report(val$df, "GSE17536_validation")

# Next steps to consider (not run): add stage / age / MMR status as covariates
# in the Cox models (columns exist in GSE39582 pData), and test RFS as well as OS.

## ---- 8. Clinical adjustment and filtered signature -------------------------
# Does the risk score add prognostic information beyond stage, age, sex (and MMR
# in GSE39582)? Stage is collapsed to I-II / III / IV in both cohorts.
stage3 <- function(x) {
  x <- trimws(as.character(x)); out <- rep(NA_character_, length(x))
  out[x %in% c("0", "1", "2")] <- "I-II"; out[x == "3"] <- "III"; out[x == "4"] <- "IV"
  factor(out, levels = c("I-II", "III", "IV"))
}
sex_f <- function(x) factor(tolower(trimws(as.character(x))), levels = c("female", "male"))

pd39 <- pData(es39582); pd17 <- pData(es17536)
mmr_raw <- trimws(as.character(pd39[["mmr.status:ch1"]]))
mmr_raw[mmr_raw %in% c("N/A", "NA", "")] <- NA
clin39582 <- data.frame(
  sample = colnames(es39582),
  stage  = stage3(pd39[["tnm.stage:ch1"]]),
  age    = suppressWarnings(as.numeric(as.character(pd39[["age.at.diagnosis (year):ch1"]]))),
  sex    = sex_f(pd39[["Sex:ch1"]]),
  mmr    = factor(mmr_raw, levels = c("pMMR", "dMMR")),
  stringsAsFactors = FALSE)
clin17536 <- data.frame(
  sample = colnames(es17536),
  stage  = stage3(pd17[["ajcc_stage:ch1"]]),
  age    = suppressWarnings(as.numeric(as.character(pd17[["age:ch1"]]))),
  sex    = sex_f(pd17[["gender:ch1"]]),
  stringsAsFactors = FALSE)

hr_row <- function(fit, term) {
  s <- summary(fit)
  sprintf("HR %.2f [%.2f-%.2f], p = %.3g",
          s$conf.int[term, "exp(coef)"], s$conf.int[term, "lower .95"],
          s$conf.int[term, "upper .95"], s$coefficients[term, "Pr(>|z|)"])
}

report_cox <- function(d, label, extra = character()) {
  vars <- c("time", "event", "grp", "risk_sd", "stage", "age", "sex", extra)
  d <- d[complete.cases(d[, vars]), ]
  rhs <- paste(c("stage", "age", "sex", extra), collapse = " + ")
  m_un <- coxph(Surv(time, event) ~ grp, data = d)
  m_cl <- coxph(as.formula(paste("Surv(time, event) ~", rhs)), data = d)
  m_ad <- coxph(as.formula(paste("Surv(time, event) ~ grp +", rhs)), data = d)
  m_ct <- coxph(as.formula(paste("Surv(time, event) ~ risk_sd +", rhs)), data = d)
  lrt  <- anova(m_cl, m_ad)
  cat(sprintf("\n%s: n = %d, events = %d (adjusted for %s)\n",
              label, nrow(d), sum(d$event), gsub(" \\+ ", ", ", rhs)))
  cat("  unadjusted  high vs low :", hr_row(m_un, "grphigh"), "\n")
  cat("  adjusted    high vs low :", hr_row(m_ad, "grphigh"), "\n")
  cat("  adjusted    per 1 SD    :", hr_row(m_ct, "risk_sd"), "\n")
  cat(sprintf("  added value of signature (LRT) p = %.3g | C-index: clinical %.3f -> +signature %.3f\n",
              lrt[2, "Pr(>|Chi|)"], summary(m_cl)$concordance[1], summary(m_ad)$concordance[1]))
  cat("  risk group by stage:\n"); print(table(risk = d$grp, stage = d$stage))
  invisible(NULL)
}

evaluate_signature <- function(cv, s, gene_set, label) {
  rk <- function(X) as.numeric(predict(cv, X[, gene_set], s = s, type = "link"))
  tr <- trn$df[, c("sample", "time", "event")]; va <- val$df[, c("sample", "time", "event")]
  tr$risk <- rk(trn$X); va$risk <- rk(val$X)
  cut <- median(tr$risk); sd_tr <- sd(tr$risk)        # fixed on training data
  for (nm in c("tr", "va")) {
    d <- get(nm)
    d$grp <- factor(ifelse(d$risk > cut, "high", "low"), levels = c("low", "high"))
    d$risk_sd <- d$risk / sd_tr
    assign(nm, d)
  }
  add_clin <- function(d, clin) {
    m <- match(d$sample, clin$sample); cbind(d, clin[m, setdiff(names(clin), "sample"), drop = FALSE])
  }
  tr <- add_clin(tr, clin39582); va <- add_clin(va, clin17536)
  cat("\n=====", label, "(", length(gene_set), "candidate genes in model ) =====")
  report_cox(va, "GSE17536 VALIDATION")
  report_cox(tr, "GSE39582 training (optimistic: model was fit here)", extra = "mmr")
  invisible(list(train = tr, valid = va))
}

res_orig <- evaluate_signature(cvcox, s_cox, uni_sig$gene, "ORIGINAL signature")

# Filtered signature: drop sex-linked genes and unannotated/antisense entries,
# then refit the penalized Cox model on the remaining univariate hits.
sex_genes <- c("KDM5D", "DDX3Y", "UTY", "EIF1AY", "RPS4Y1", "RPS4Y2", "ZFY", "USP9Y",
               "TXLNGY", "NLGN4Y", "TMSB4Y", "XIST", "TSIX", "PRKY", "TBL1Y", "AMELY")
unannot <- grepl("^LOC[0-9]|^LINC|-AS[0-9]+$|^MIR[0-9]|^SNOR", uni_sig$gene)
cand2 <- setdiff(uni_sig$gene[!unannot], sex_genes)
cat("\nRemoved", length(uni_sig$gene) - length(cand2), "genes:",
    paste(setdiff(uni_sig$gene, cand2), collapse = ", "), "\n")

set.seed(2026)
cv2 <- cv.glmnet(trn$X[, cand2], ysurv, family = "cox", alpha = 0.9, nfolds = 10)
s2 <- "lambda.1se"; b2 <- as.matrix(coef(cv2, s = s2))
if (sum(b2[, 1] != 0) == 0) { message("filtered model: lambda.1se empty; using lambda.min")
  s2 <- "lambda.min"; b2 <- as.matrix(coef(cv2, s = s2)) }
genes2 <- rownames(b2)[b2[, 1] != 0]
write.csv(data.frame(gene = genes2, coef = b2[genes2, 1]),
          "results/prognostic_signature_filtered.csv", row.names = FALSE)
res_filt <- evaluate_signature(cv2, s2, cand2, paste0("FILTERED signature (", s2, ", ", length(genes2), " genes)"))

## ---- 9. Freeze the primary model (do this BEFORE touching TCGA) -----------
# Primary model = FILTERED signature. After this section nothing is re-fit,
# re-tuned or re-cut; TCGA is only scored with what is saved here.
dir.create("results/frozen_model", showWarnings = FALSE)

beta <- b2[genes2, 1]                                   # non-zero coefficients only
manual_score <- function(X) as.numeric(X[, genes2, drop = FALSE] %*% beta)
tr_score <- manual_score(trn$X)
va_score <- manual_score(val$X)

# Sanity check 1: our hand-computed score matches glmnet's own linear predictor
chk <- cor(tr_score, as.numeric(predict(cv2, trn$X[, cand2], s = s2, type = "link")))
stopifnot(chk > 0.999)

# Sanity check 2: the frozen score reproduces the validation association
m_chk <- coxph(Surv(val$df$time, val$df$event) ~ I(va_score / sd(tr_score)))
cat("\nFrozen score, GSE17536 unadjusted HR per training-SD:",
    sprintf("%.2f [%.2f-%.2f]\n", summary(m_chk)$conf.int[1, 1],
            summary(m_chk)$conf.int[1, 3], summary(m_chk)$conf.int[1, 4]))

frozen <- list(
  created   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  model     = "Penalized Cox (glmnet, alpha = 0.9), filtered gene set",
  lambda    = s2, genes = genes2, coef = beta,
  training  = list(cohort = "GSE39582 tumors with OS", n = nrow(trn$df),
                   events = sum(trn$df$event)),
  train_score_median = median(tr_score), train_score_sd = sd(tr_score),
  scoring   = paste("log2 expression; z-score EACH GENE within the cohort being scored;",
                    "score = sum(coef * z)"),
  removed_genes = c(sex_genes, grep("^LOC[0-9]|^LINC|-AS[0-9]+$|^MIR[0-9]|^SNOR",
                                    uni_sig$gene, value = TRUE)),
  validation_gse17536 = list(n = nrow(val$df), events = sum(val$df$event))
)
saveRDS(frozen, "results/frozen_model/frozen_signature.rds")
write.csv(data.frame(gene = genes2, coef = beta),
          "results/frozen_model/signature_coefficients.csv", row.names = FALSE)

plan <- c(
  "PRE-SPECIFIED VALIDATION PLAN (written BEFORE any TCGA outcome data were analysed)",
  paste("Date frozen:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "MODEL: frozen_signature.rds / signature_coefficients.csv (no re-fitting, no re-selection).",
  "COHORT: TCGA-COAD primary tumours. Exclude rectal (READ) cases, missing OS, OS time <= 0.",
  "OUTCOME: overall survival (use the curated TCGA-CDR clinical resource).",
  "SCORING: log2(normalized expression + 1); z-score each gene within the TCGA tumours;",
  "         score = sum(coefficient x z). Report how many of the signature genes are present;",
  "         if fewer than 90% are present, STOP and discuss before proceeding.",
  "",
  "PRIMARY ANALYSIS: Cox model, hazard ratio per 1 SD of the TCGA score,",
  "         adjusted for stage (I-II / III / IV), age at diagnosis, sex.",
  "PRIMARY SUCCESS CRITERION: [CONFIRM WITH ADVISORS BEFORE RUNNING]",
  "         default = adjusted HR > 1 with the 95% CI excluding 1.",
  "",
  "SECONDARY: high vs low by TCGA-cohort median; added value over stage+age+sex",
  "         (likelihood-ratio test, C-index); MSI/MMR-adjusted model if MSI status is available.",
  "SENSITIVITY: follow-up capped at 60 months; excluding stage IV.",
  "",
  "NOT ALLOWED after seeing TCGA outcomes: changing genes, coefficients, cutoffs, covariates,",
  "or the success criterion. Any deviation is reported as post-hoc."
)
writeLines(plan, "results/frozen_model/TCGA_validation_plan.txt")

fz <- list.files("results/frozen_model", full.names = TRUE)
fz <- fz[basename(fz) != "checksums.txt"]
writeLines(paste(tools::md5sum(fz), basename(fz)), "results/frozen_model/checksums.txt")
cat("\nFrozen model saved to results/frozen_model (", length(genes2), "genes ).\n",
    "Share/commit this folder BEFORE running TCGA so the timestamp proves the order.\n")

## ---- Session info ----------------------------------------------------------
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat("\nDone. Outputs are in ./results\n")