# Plot sequencing depth to make sure I removed all samples with low # of reads
seq_depth_16S<- data.table(as(sample_data(physeq_16S), "data.frame"),
                       TotalReads = sample_sums(physeq_16S), keep.rownames = TRUE)
setnames(seq_depth_16S, "rn", "SampleID")
plot_seq_depth_16S <- ggplot(seq_depth_16S, aes(TotalReads)) + geom_histogram() + 
  ggtitle("Sequencing Depth 16S")
plot_seq_depth_16S