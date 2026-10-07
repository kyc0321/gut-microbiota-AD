## _helpers.R — shared utilities
## JWT handling: prefer an already-set OPENGWAS_JWT (e.g. from ~/.Renviron,
## which is kept current) over the local cache file. The cache file path is
## resolved relative to the working directory (project root) rather than a
## hardcoded absolute path, since the project folder has been relocated at
## least once (was previously directly under Desktop, now nested under the
## "#3 periodontitis..." folder) and a hardcoded path silently breaks then.
JWT_FILE <- file.path("data", ".jwt")
if (nzchar(Sys.getenv("OPENGWAS_JWT"))) {
  # Already set (e.g. via ~/.Renviron) -- use it as-is, don't override.
} else if (file.exists(JWT_FILE)) {
  Sys.setenv(OPENGWAS_JWT = readLines(JWT_FILE, warn = FALSE)[1])
} else {
  stop("No OPENGWAS_JWT found in environment and no cache file at ", JWT_FILE,
       ". Set OPENGWAS_JWT in ~/.Renviron or place a valid token at this path.")
}

delta_se <- function(a, sa, b, sb) sqrt(a^2 * sb^2 + b^2 * sa^2)
