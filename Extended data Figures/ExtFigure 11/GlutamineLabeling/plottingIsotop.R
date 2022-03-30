plottingIsotop <- function(dataframeIsotops, columns){

	tbl_gln <- cbind(dataframeIsotops[, ..columns], cell_line=dataframeIsotops$cell_line)
	tbl_gln[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=columns]
	avrg_sum <- mean(tbl_gln$sample_sum, na.rm = TRUE)

	tbl_gln_norm <- cbind((tbl_gln[, ..columns]/tbl_gln$sample_sum)*(tbl_gln$sample_sum/avrg_sum), cell_line=dataframeIsotops$cell_line)
	tbl_gln_norm_melt <- melt.data.table(tbl_gln_norm)
  	tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.infinite(tbl_gln_norm_melt$value), ]
  	tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.na(tbl_gln_norm_melt$value), ]
	tbl_gln_norm_melt[, isotope := "M+0"]
	tbl_gln_norm_melt[grep("C1", tbl_gln_norm_melt$variable), isotope := "M+1"]
	tbl_gln_norm_melt[grep("C2", tbl_gln_norm_melt$variable), isotope := "M+2"]
	tbl_gln_norm_melt[grep("C3", tbl_gln_norm_melt$variable), isotope := "M+3"]
	tbl_gln_norm_melt[grep("C4", tbl_gln_norm_melt$variable), isotope := "M+4"]
	tbl_gln_norm_melt[grep("C5", tbl_gln_norm_melt$variable), isotope := "M+5"]
	tbl_gln_norm_melt[grep("C6", tbl_gln_norm_melt$variable), isotope := "M+6"]

	df.summary <- tbl_gln_norm_melt %>%
	  group_by(cell_line, isotope) %>%
	  summarise(
	    sd = sd(value, na.rm = TRUE),
	    len = mean(value, na.rm = TRUE)
	  )
	df.summary <- df.summary %>%
	  mutate(pos = cumsum(len))

	gln_plt_poolsize <- 
		ggplot(tbl_gln_norm_melt, aes(x=cell_line, y=value, fill=(isotope))) +
		  geom_bar(position = position_stack(reverse = TRUE), stat="summary", fun="mean", alpha=0.6) +
		  geom_errorbar(data=df.summary, aes(x=cell_line, ymin=pos, ymax=pos+sd, color=(isotope)), inherit.aes = FALSE, width=0.3, size = 0.5) +
		  labs(title=paste0(columns[1]), fill="Isotopologues", color="Isotopologues", y="Normalized incorporation") +
		  scale_color_npg()+
		  scale_fill_npg()+
		  theme_test() +
		  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color = "black"), axis.text.y=element_text(color = "black"), axis.title.y=element_text(size=10), plot.title=element_text(size=10), legend.position = "bottom", legend.title=element_blank())

	return(gln_plt_poolsize)

}


plottingIsotopGrouped <- function(dataframeIsotops, columns){

	tbl_gln <- cbind(dataframeIsotops[, ..columns], cell_name=dataframeIsotops$cell_name)
	tbl_gln[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=columns]
	avrg_sum <- mean(tbl_gln$sample_sum, na.rm = TRUE)

	tbl_gln_norm <- cbind((tbl_gln[, ..columns]/tbl_gln$sample_sum)*(tbl_gln$sample_sum/avrg_sum), cell_name=dataframeIsotops$cell_name)
	tbl_gln_norm_melt <- melt.data.table(tbl_gln_norm)
  	tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.infinite(tbl_gln_norm_melt$value), ]
  	tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.na(tbl_gln_norm_melt$value), ]
	tbl_gln_norm_melt[, isotope := "M+0"]
	tbl_gln_norm_melt[grep("C1", tbl_gln_norm_melt$variable), isotope := "M+1"]
	tbl_gln_norm_melt[grep("C2", tbl_gln_norm_melt$variable), isotope := "M+2"]
	tbl_gln_norm_melt[grep("C3", tbl_gln_norm_melt$variable), isotope := "M+3"]
	tbl_gln_norm_melt[grep("C4", tbl_gln_norm_melt$variable), isotope := "M+4"]
	tbl_gln_norm_melt[grep("C5", tbl_gln_norm_melt$variable), isotope := "M+5"]
	tbl_gln_norm_melt[grep("C6", tbl_gln_norm_melt$variable), isotope := "M+6"]

	df.summary <- tbl_gln_norm_melt %>%
	  group_by(cell_name, isotope) %>%
	  summarise(
	    sd = sd(value, na.rm = TRUE),
	    len = mean(value, na.rm = TRUE)
	  )
	df.summary <- df.summary %>%
	  mutate(pos = cumsum(len))

	gln_plt_poolsize <- 
		ggplot(tbl_gln_norm_melt, aes(x=cell_name, y=value, fill=(isotope))) +
		  geom_bar(position = position_stack(reverse = TRUE), stat="summary", fun="mean", alpha=0.6) +
		  geom_errorbar(data=df.summary, aes(x=cell_name, ymin=pos, ymax=pos+sd, color=(isotope)), inherit.aes = FALSE, width=0.3, size = 0.5) +
		  labs(title=paste0(columns[1]), fill="Isotopologues", color="Isotopologues", y="Normalized incorporation") +
		  scale_color_npg()+
		  scale_fill_npg()+
		  theme_test() +
		  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color = "black"), axis.text.y=element_text(color = "black"), axis.title.y=element_text(size=10), plot.title=element_text(size=10), legend.position = "bottom", legend.title=element_blank())

	  return(gln_plt_poolsize)

}





plottingTICgrouped <- function(dataframeIsotops, columns){

  tbl_gln <- cbind(dataframeIsotops[, ..columns], cell_name=dataframeIsotops$cell_name, protein=dataframeIsotops$protein)
  tbl_gln[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=columns]
  tbl_gln[, ioncountprotein := sample_sum/protein]
  tbl_gln_norm_melt <- melt.data.table(tbl_gln[, c("cell_name", "ioncountprotein")])
  tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.infinite(tbl_gln_norm_melt$value), ]
  tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.na(tbl_gln_norm_melt$value), ]

  df.summary <- tbl_gln_norm_melt %>%
    group_by(cell_name) %>%
    summarise(
      sd = sd(value, na.rm = TRUE),
      len = mean(value, na.rm = TRUE)
    )

  gln_plt_poolsize <- 
    ggplot(tbl_gln_norm_melt, aes(x=cell_name, y=value, fill=cell_name)) +
    geom_errorbar(data=df.summary, aes(x=cell_name, ymin=len-sd, ymax=len+sd), inherit.aes = FALSE, width=0.3, size = 0.5) +
    geom_bar(stat="summary", fun="mean",) +
    geom_jitter(shape=21, width = 0.2, fill="white", color="black") +
    labs(title=paste0(columns[1]),y="Total ion count / protein amount") +
    scale_fill_aaas() +
    theme_test() +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color = "black"), axis.text.y=element_text(color = "black"), axis.title.y=element_text(size=10), plot.title=element_text(size=10), legend.position = "right", legend.title=element_blank())

	return(gln_plt_poolsize)

}









plottingFractionsGrouped <- function(dataframeIsotops, columns){

	tbl_gln <- cbind(dataframeIsotops[, ..columns], cell_name=dataframeIsotops$cell_name, protein=dataframeIsotops$protein)
	tbl_gln[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=columns]
	avrg_sum <- mean(tbl_gln$sample_sum, na.rm = TRUE)

	tbl_gln_norm <- cbind((tbl_gln[, ..columns]/tbl_gln$sample_sum), cell_name=dataframeIsotops$cell_name) # could additionally normalize to overall signal: *(tbl_gln$sample_sum/avrg_sum)
	tbl_gln_norm_melt <- melt.data.table(tbl_gln_norm)
	tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.infinite(tbl_gln_norm_melt$value), ]
	tbl_gln_norm_melt <- tbl_gln_norm_melt[!is.na(tbl_gln_norm_melt$value), ]
	tbl_gln_norm_melt[, isotope := "M+0"]
	tbl_gln_norm_melt[grep("C1", tbl_gln_norm_melt$variable), isotope := "M+1"]
	tbl_gln_norm_melt[grep("C2", tbl_gln_norm_melt$variable), isotope := "M+2"]
	tbl_gln_norm_melt[grep("C3", tbl_gln_norm_melt$variable), isotope := "M+3"]
	tbl_gln_norm_melt[grep("C4", tbl_gln_norm_melt$variable), isotope := "M+4"]
	tbl_gln_norm_melt[grep("C5", tbl_gln_norm_melt$variable), isotope := "M+5"]
	tbl_gln_norm_melt[grep("C6", tbl_gln_norm_melt$variable), isotope := "M+6"]

	df.summary <- tbl_gln_norm_melt %>%
	group_by(cell_name, isotope) %>%
	summarise(
	sd = sd(value, na.rm = TRUE),
	len = mean(value, na.rm = TRUE)
	)
	df.summary <- df.summary %>%
	group_by(cell_name) %>%
	mutate(
	sum = sum(len, na.rm = TRUE)
	)
	df.summary$norm_len <- df.summary$len / df.summary$sum

	library(RColorBrewer)
	colfunc <- colorRampPalette(brewer.pal(9,"BuPu"))(9)

	gln_plt_poolsize <- 
	ggplot(df.summary, aes(x=cell_name, y=norm_len, fill=(isotope)))+
	geom_col(position = position_stack(reverse = TRUE), alpha = 0.8)+
	labs(title=paste0(columns[1]),y="Relative fractional incorporation")+
	scale_fill_manual(values=colfunc[3:9], guide = guide_legend(reverse = TRUE))+
	theme_test()+
	theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color = "black"), axis.text.y=element_text(color = "black"), axis.title.y=element_text(size=10), plot.title=element_text(size=10), legend.position = "right", legend.title=element_blank())


	return(gln_plt_poolsize)

}




plottingFractionsRatio <- function(dataframeIsotops, columns, numeratorIsotop, denominatorIsotop){

	tbl_gln <- cbind(dataframeIsotops[, ..columns], cell_name=dataframeIsotops$cell_name, protein=dataframeIsotops$protein)
	tbl_gln[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=columns]
	avrg_sum <- mean(tbl_gln$sample_sum, na.rm = TRUE)

	tbl_gln_norm <- cbind((tbl_gln[, ..columns]/tbl_gln$sample_sum), cell_name=dataframeIsotops$cell_name) # could additionally normalize to overall signal: *(tbl_gln$sample_sum/avrg_sum)
	
	ifelse(paste0(numeratorIsotop) == "C0" | paste0(denominatorIsotop) == "C0",
		ifelse(paste0(numeratorIsotop) == "C0",
			tbl_gln_norm[ , ratio := (tbl_gln_norm[, 1, with=FALSE] / tbl_gln_norm[, grep(paste0(denominatorIsotop), colnames(tbl_gln_norm)), with=FALSE])],
			tbl_gln_norm[ , ratio := (tbl_gln_norm[, grep(paste0(numeratorIsotop), colnames(tbl_gln_norm)), with=FALSE] / tbl_gln_norm[, 1, with=FALSE])]
			),
		tbl_gln_norm[ , ratio := (tbl_gln_norm[, grep(paste0(numeratorIsotop), colnames(tbl_gln_norm)), with=FALSE] / tbl_gln_norm[, grep(paste0(denominatorIsotop), colnames(tbl_gln_norm)), with=FALSE])]
		)

	gln_plt_poolsize <- 
	ggplot(tbl_gln_norm, aes(x=cell_name, y=ratio, fill=cell_name))+
	  geom_boxplot(outlier.shape = NA, alpha = 0.8)+
	  geom_jitter(shape=21, width = 0.2, fill="white", color="black") +
	  labs(title=paste0(columns[1]), y=paste0(numeratorIsotop,"/",denominatorIsotop))+
	  scale_fill_aaas() +
	  theme_test()+
	  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color = "black"), axis.text.y=element_text(color = "black"), axis.title.y=element_text(size=10), plot.title=element_text(size=10), legend.position = "none", legend.title=element_blank())

	return(gln_plt_poolsize)

}





plottingMetaboliteLabellingBox <- function(dataframeIsotops, columns){

	tbl_gln <- cbind(dataframeIsotops[, ..columns], cell_name=dataframeIsotops$cell_name, protein=dataframeIsotops$protein)
	tbl_gln[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=columns]
	avrg_sum <- mean(tbl_gln$sample_sum, na.rm = TRUE)
	number_isotop <- length(columns)

	label_list <- list()
	for (i in 1:number_isotop) {
	  label_list[[i]] <- tbl_gln[, i, with=FALSE]*((i-1)/(number_isotop-1))
	}
	tbl_gln2 <- data.table(do.call(cbind, label_list), sample_sum=tbl_gln$sample_sum)

	tbl_gln_norm <- cbind((tbl_gln2[, ..columns]/tbl_gln$sample_sum), cell_name=tbl_gln$cell_name)
	tbl_gln_norm[, labelling := rowSums(.SD, na.rm = FALSE), .SDcols=columns]

	gln_plt_poolsize <- 
	ggplot(tbl_gln_norm, aes(x=cell_name, y=labelling, fill=cell_name))+
	  geom_boxplot(outlier.shape = NA, alpha = 0.8)+
	  geom_jitter(shape=21, width = 0.2, fill="white", color="black") +
	  labs(title=paste0(columns[1]), y="Normalized labelling")+
	  scale_fill_aaas() +
	  theme_test()+
	  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color = "black"), axis.text.y=element_text(color = "black"), axis.title.y=element_text(size=10), plot.title=element_text(size=10), legend.position = "none", legend.title=element_blank())

	return(gln_plt_poolsize)

}





plottingMetaboliteLabelling <- function(dataframeIsotops, columns){

	tbl_gln <- cbind(dataframeIsotops[, ..columns], cell_name=dataframeIsotops$cell_name, protein=dataframeIsotops$protein)
	tbl_gln[, sample_sum := rowSums(.SD, na.rm = FALSE), .SDcols=columns]
	avrg_sum <- mean(tbl_gln$sample_sum, na.rm = TRUE)
	number_isotop <- length(columns)

	label_list <- list()
	for (i in 1:number_isotop) {
	  label_list[[i]] <- tbl_gln[, i, with=FALSE]*((i-1)/(number_isotop-1))
	}
	tbl_gln2 <- data.table(do.call(cbind, label_list), sample_sum=tbl_gln$sample_sum)

	tbl_gln_norm <- cbind((tbl_gln2[, ..columns]/tbl_gln$sample_sum), cell_name=tbl_gln$cell_name)
	tbl_gln_norm[, labelling := rowSums(.SD, na.rm = FALSE), .SDcols=columns]

	tbl_gln_norm_sel <- tbl_gln_norm[, c("cell_name", "labelling")]
	# tbl_gln_norm_sel <- tbl_gln_norm_sel[!is.infinite(tbl_gln_norm_sel$value), ]
	# tbl_gln_norm_sel <- tbl_gln_norm_sel[!is.na(tbl_gln_norm_sel$value), ]

	df.summary <- tbl_gln_norm_sel %>%
	group_by(cell_name) %>%
	summarise(
	  sd = sd(labelling, na.rm = TRUE),
	  len = mean(labelling, na.rm = TRUE)
	)

	gln_plt_poolsize <- 
	ggplot(tbl_gln_norm, aes(x=cell_name, y=labelling, fill=cell_name))+
		geom_errorbar(data=df.summary, aes(x=cell_name, ymin=len-sd, ymax=len+sd), inherit.aes = FALSE, width=0.3, size = 0.5) +
		geom_bar(stat="summary", fun="mean",) +
		geom_jitter(shape=21, width = 0.2, fill="white", color="black") +
		labs(title=paste0(columns[1]), y="Normalized labelling")+
		scale_fill_aaas() +
		theme_test()+
		theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color = "black"), axis.text.y=element_text(color = "black"), axis.title.y=element_text(size=10), plot.title=element_text(size=10), legend.position = "none", legend.title=element_blank())

	return(gln_plt_poolsize)

}



