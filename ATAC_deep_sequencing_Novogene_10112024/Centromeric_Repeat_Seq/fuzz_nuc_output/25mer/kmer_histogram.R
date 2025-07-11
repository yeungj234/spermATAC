setwd("/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/count_kmers_in_fastq_verbose/output")
kmer_summary <- dir(pattern="summary.tsv")
total_reads <- c(66233269, 91756371, 81105822 , 54252756, 62717380, 115703381)
names(total_reads) <- gsub("_counts_summary.tsv", "", kmer_summary)
hist <- list()
for(i in 1:length(kmer_summary)){
        hist[[i]] <- read.table(kmer_summary[i], header=T, stringsAsFactors = FALSE)
        hist[[i]]<- hist[[i]][!hist[[i]]$kmer_frequencies %in% "num_kmers_matched", ]
        hist[[i]]$kmer_frequencies <- as.numeric(hist[[i]]$kmer_frequencies)
        hist[[i]]$Percent_total_reads <- (hist[[i]]$Number_of_reads/total_reads[i])*100
        hist[[i]]$sample <- names(total_reads)[i]
}
names(hist) <- names(total_reads) 
hist <- hist[c(1,3,4,2,5,6)]

hist_merged <- do.call(rbind, hist)
hist_merged$sample <- factor(hist_merged$sample, levels=names(hist))
hist_merged$Condition <- ifelse(hist_merged$sample %in% names(hist)[1:3], "Interphase", "Mitosis")
library(ggplot2)
pdf_dir <- "/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Centromeric_Repeat_Seq/count_kmers_in_fastq_verbose/Figures"
pdf(paste0(pdf_dir, "/IvsM_25mer_freq_in_dedup_fastq-reads.pdf"), width=8, height=4)
ggplot(hist_merged[!hist_merged$kmer_frequencies == 0, ], aes(x=kmer_frequencies, y=Percent_total_reads, fill=sample, alpha=0.3), color="black")+
        geom_bar(stat="identity", position="identity",  color="black", size=0.1)+
        scale_fill_manual(values=c("#BADEFAFF","#90CAF8FF","#64B4F6FF", "#9FA7D9FF", "#5B6BBFFF", "#19227EFF"))+
        ylim(0,0.03)+
        facet_wrap(~Condition)+
        ylab("% Total Reads")+
        xlab("Number of kmers found")+
        theme_minimal()

ggplot(hist_merged[!hist_merged$kmer_frequencies == 0, ], aes(x=kmer_frequencies, y=Percent_total_reads, fill=sample))+
        geom_bar(stat="identity", size=0.1,  color="black")+
        scale_fill_manual(values=c("#BADEFAFF","#90CAF8FF","#64B4F6FF", "#9FA7D9FF", "#5B6BBFFF", "#19227EFF"))+
        ylim(0,0.03)+
        facet_wrap(Condition~sample)+
        ylab("% Total Reads")+
        xlab("Number of kmers found")+
        theme_minimal()

ggplot(hist_merged[!hist_merged$kmer_frequencies == 0, ], aes(x=kmer_frequencies, y=Percent_total_reads, fill=sample))+
        geom_bar(stat="identity", size=0.05, color="black")+
        scale_fill_manual(values=c("#BADEFAFF","#90CAF8FF","#64B4F6FF", "#9FA7D9FF", "#5B6BBFFF", "#19227EFF"))+
        ylim(0,0.03)+
        facet_wrap(Condition~sample)+
        ylab("% Total Reads")+
        xlab("Number of kmers found")+
        theme_minimal()

ggplot(hist_merged[!hist_merged$kmer_frequencies == 0, ], aes(x=kmer_frequencies, y=Percent_total_reads, color=sample,))+
        geom_line(stat="identity", position="identity", size=1)+
        scale_color_manual(values=c("#BADEFAFF","#90CAF8FF","#64B4F6FF", "#9FA7D9FF", "#5B6BBFFF", "#19227EFF"))+
        ylim(0,0.03)+
        ylab("% Total Reads")+
        xlab("Number of kmers found")+
        theme_minimal()
dev.off()

Reads_with_kmers_stats <- lapply(hist, function(x){percent_no_kmers_reads <- x[x$kmer_frequencies ==0, ]$Percent_total_reads
                        percent_kmer_reads <- 100-percent_no_kmers_reads
                        num_kmer_reads <- sum(x[!x$kmer_frequencies ==0, ]$Number_of_reads)
                        df <- data.frame(Percent_kmer_reads=percent_kmer_reads, Number_of_kmer_reads=num_kmer_reads)
                        return(df)
                        })
Reads_with_kmers_stats <- do.call(rbind, Reads_with_kmers_stats)
Reads_with_kmers_stats$total_reads <- total_reads[c(1,3,4,2,5,6)]
