# Wilcoxon test
processFile = function(filepath) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    new_line = unlist(strsplit(line, split = "\t"))
    if ( length(line) == 0 ) {
      break
    }
    copro = as.numeric(new_line[1:8])
    soil = as.numeric(new_line[793:795])
    traditional = as.numeric(new_line[9:376])
    urban = as.numeric(new_line[377:792])
    modern = as.numeric(new_line[9:792])
    
    res1 = wilcox.test(copro, traditional, alternative=c("two.sided"))
    res2 = wilcox.test(copro, urban, alternative=c("two.sided"))
    res3 = wilcox.test(traditional, urban, alternative=c("two.sided"))
    res4 = wilcox.test(soil, modern, alternative=c("two.sided"))
    
    write.table(res1$p.value, file = "/Users/marshacw/Documents/wilcoxon/res1.csv", append = TRUE, sep = "\t")
    write.table(res2$p.value, file = "/Users/marshacw/Documents/wilcoxon/res2.csv", append = TRUE, sep = "\t")
    write.table(res3$p.value, file = "/Users/marshacw/Documents/wilcoxon/res3.csv", append = TRUE, sep = "\t")
    write.table(res4$p.value, file = "/Users/marshacw/Documents/wilcoxon/res4.csv", append = TRUE, sep = "\t")
  }
  close(con)}

processFile("/Users/marshacw/Documents/for_wilcoxon.txt")

# FDR correction for p-values
phylum_wilcoxon= read.table("/Users/marshacw/Documents/wilcoxon_results.txt",sep="\t",header=T, stringsAsFactors = F, quote="", row.names=1)

p_adjust_copro_vs_trad<-as.data.frame(p.adjust(phylum_wilcoxon[,1], method = "fdr"))
p_adjust_copro_vs_urban<-as.data.frame(p.adjust(phylum_wilcoxon[,2], method = "fdr"))
p_adjust_trad_vs_urban<-as.data.frame(p.adjust(phylum_wilcoxon[,3], method = "fdr"))
p_adjust_soil_vs_modern<-as.data.frame(p.adjust(phylum_wilcoxon[,4], method = "fdr"))

final_p_adjust <- paste(row.names(phylum_wilcoxon), "\t",p_adjust_copro_vs_trad[,1], "\t",p_adjust_copro_vs_urban[,1], "\t",p_adjust_trad_vs_urban[,1], "\t",p_adjust_soil_vs_modern[,1])

write.table(final_p_adjust, "/Users/marshacw/Documents/wilcoxon_padj.txt", sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE, quote=FALSE)
