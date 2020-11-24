# Fisher test
fish <- read.table("/Users/marshacw/Documents/for_fisher.txt", header = TRUE, sep = "\t", row.names=1)
result = apply(as.matrix(fish[,1:4]), 1, function(x)
  fisher.test(matrix(x, ncol = 2), alternative = "two.sided", conf.int=TRUE)$p.value)

write.table(result, "/Users/marshacw/Documents/fisher_results.txt", sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE, quote=FALSE)