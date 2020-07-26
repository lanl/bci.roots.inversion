

rl <- readLines("manuscript.rmd")

library(stringr)
ref <- str_extract_all(rl, "@[a-zA-Z0-9]++")
length(unique(do.call(c, ref)))
