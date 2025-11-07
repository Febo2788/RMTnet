user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
.libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({ library(ALL); library(Biobase) })
data(ALL)
mat_full <- exprs(ALL)
iqr <- apply(mat_full, 1, IQR)
probes <- rownames(mat_full)[order(iqr, decreasing = TRUE)[1:1000]]
write.csv(data.frame(probe_id = probes), "C:/Users/felix/downloads/test/probes.csv", row.names = FALSE)
