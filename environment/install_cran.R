# Packages missing from conda-forge for osx-arm64; built from CRAN source inside vg-local.
pkgs <- c("phylolm", "corHMM", "geiger", "castor")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org", type = "source", Ncpus = 8)
for (p in pkgs) cat(p, as.character(packageVersion(p)), "\n")
