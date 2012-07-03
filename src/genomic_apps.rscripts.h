const char *RSCRIPT_TEMPLATE_HEATMAP = \
"\n##\n## USAGE: genomic_apps.heatmap.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\nnorm_rows <- function(X) { for (i in 1:nrow(X)) X[i,] <- (X[i,]-mean(X[i,]))/sd(X[i,]); return(X); }\n\nargs <- commandArgs(trailingOnly=T);\ndata_file <- args[1];\nparam_file <- args[2];\nimage_file <- args[3];\n\nparams <- readLines(param_file);\nshift_upstream <- as.numeric(params[1]);\nshift_downstream <- as.numeric(params[2]);\nheatmap_colors <- strsplit(params[3],',')[[1]];\nheatmap_title <- strsplit(params[4],',')[[1]];\nheatmap_xlab <- params[5];\nheatmap_ylab <- params[6];\nheatmap_size <- as.numeric(strsplit(params[7],',')[[1]]);\nheatmap_resolution <- as.numeric(params[8]);\nn_heatmaps <- as.numeric(params[9]);\nimage_type <- rev(unlist(strsplit(image_file,'.',fixed=T)))[1];\n\n\nlibrary('MASS');\nlibrary('preprocessCore');           # from Bioconductor\nlibrary('gplots');                   # from Bioconductor\n\n# load data\nD <- as.matrix(read.table(data_file,row.names=1,sep='\t'));\n  \n# create combined heatmap (main version)\nif (image_type=='tif') {\n  tiff(image_file,width=heatmap_size[1],height=heatmap_size[2],res=heatmap_resolution,compression='lzw');\n} else if (image_type=='pdf') {\n  pdf(image_file); \n}\n\npar(fig=c(0,1,0,1),mar=c(2,2,0,0)); \nplot.new();\nmtext(heatmap_xlab,side=1);\nmtext(heatmap_ylab,side=2);\n\n\nd <- ncol(D)/n_heatmaps;\nI <- 1:d;\ndj <- 0.9/n_heatmaps;\nfor (j in 1:n_heatmaps) {\n  par(fig=c(0.1+(j-1)*dj,0.1+j*dj,0.1,1),mar=c(0.5,0.5,3,0.5),new=TRUE);\n  colorscale <- c(colorpanel(20,low='white',high='white'),colorpanel(50,low='white',high=heatmap_colors[j]),colorpanel(30,low=heatmap_colors[j],high=heatmap_colors[j]));\n  image(z=t(norm_rows(D[,I])),col=colorscale,main=heatmap_title[j],xlab=heatmap_xlab,ylab=heatmap_ylab,xaxt='n',yaxt='n');\n  I <- I+d;\n}\n\ndev.off();\n\n";

const char *RSCRIPT_TEMPLATE_PEAKDIFF = \
"# source('c:/Aris/Research/Code/genomic_tools/genomic_apps.peakdiff.r');\n\n\n##\n## USAGE: genomic_apps.peakdiff.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\n\nsample_poisson <- function(x)\n{\n  y <- x;\n  i <- 1;\n  while (i <= length(x)) { \n    y[i] <- rpois(1,x[i]); \n	i <- i + 1;\n  }\n  return(y);\n}\n\n\n\nnormalize_matrix <- function(D)\n{\n  D_norm <- normalize.quantiles(D);\n  dimnames(D_norm) <- dimnames(D);\n  return(D_norm);\n}\n\n\nmax_abs <- function(a,b) \n{\n  if (abs(a)>abs(b)) return(a) else return(b);\n}\n\nscore_binomial <- function(x,y,n) \n{\n  x <- min(x,1);\n  y <- min(y,1);\n  a <- -pbinom(n*x,n,y,lower.tail=F,log.p=T);\n  b <- pbinom(n*x,n,y,log.p=T);\n  return(max_abs(a,b));\n}\n\nscore_fold <- function(x,y) \n{\n  return(x/y);\n}\n\nscore <- function(x,y,n,logpval) \n{\n  p <- score_binomial(x,y,n);\n  if (abs(p) < -logpval) return(1.0) else return(score_fold(x,y));\n}\n\n\ncalc_fdr_cutoff <- function(pos,neg,fdr) \n{\n  if (fdr<=0) return(Inf);\n  pos <- sort(pos);\n  neg <- sort(neg);\n  kpos <- 1;\n  kneg <- 1;\n  while ((kpos<=length(pos))&(kneg<=length(neg))) {\n    if ((length(neg)-kneg+1)/(length(pos)-kpos+1)<=fdr) { break; }\n    if (pos[kpos]<neg[kneg]) { kpos <- kpos+1; }\n    else if (pos[kpos]>neg[kneg]) { kneg <- kneg+1; }\n    else { kpos <- kpos+1; kneg <- kneg+1; }\n  }\n  if (kpos>length(pos)) { y <- Inf; }\n  else { y <- pos[kpos]; }\n  return(y);\n}\n\n\ncalc_fdr_cutoff_with_bins <- function(value,pos,neg,fdr,bin_width) \n{\n  n <- length(value);\n  t_score <- seq(0,0,length.out=n);\n  I <- order(value);\n  i <- 1; \n  value_bin <- c();\n  t_cutoff <- c();\n  while (i <= n) {\n    a <- I[i:min(i+bin_width-1,n)];\n    t_bin_cutoff <- calc_fdr_cutoff(pos[a],neg[a],fdr);\n    t_score[a] <- pos[a]/t_bin_cutoff; \n    t_cutoff <- c(t_cutoff,t_bin_cutoff,t_bin_cutoff);\n	value_bin <- c(value_bin,min(value[a]),max(value[a]));\n    i <- i + bin_width;\n  }\n  return(list(value_bin=value_bin,t_cutoff=t_cutoff,t_score=t_score));\n}\n\n\ndiff_peaks.calc <- function(D,signal_cols,ref_cols,fdr,n_fdr_bins,win_size,logpval)\n{\n  cat('Computing scores...',fill=T);\n  t_over <- apply(D,1,function(x) score(mean(x[signal_cols]),mean(x[ref_cols]),win_size,logpval));\n  t_under <- apply(D,1,function(x) score(mean(x[ref_cols]),mean(x[signal_cols]),win_size,logpval));\n  cat('Computing randomized scores...',fill=T);\n  t_sig_control <- apply(D,1,function(x) mean(c(score(x[signal_cols[1]],x[signal_cols[2]],win_size,logpval),score(x[signal_cols[2]],x[signal_cols[1]],win_size,logpval))));\n  t_ref_control <- apply(D,1,function(x) mean(c(score(x[ref_cols[1]],x[ref_cols[2]],win_size,logpval),score(x[ref_cols[2]],x[ref_cols[1]],win_size,logpval))));\n  cat('Computing FDR scores...',fill=T);\n  bin_width <- round(nrow(D)/n_fdr_bins);\n  t_over_cutoff <- calc_fdr_cutoff_with_bins(apply(D[,ref_cols],1,mean),t_over,t_ref_control,fdr,bin_width);\n  t_under_cutoff <- calc_fdr_cutoff_with_bins(apply(D[,signal_cols],1,mean),t_under,t_sig_control,fdr,bin_width);\n  cat('significant entries (over) = '); cat(sum(t_over_cutoff$t_score>=1),fill=T);\n  cat('significant entries (under) = '); cat(sum(t_under_cutoff$t_score>=1),fill=T);\n  peaks <- list(D=D,signal_cols=signal_cols,ref_cols=ref_cols,fdr=fdr,t_over=t_over,t_under=t_under,t_over_cutoff=t_over_cutoff,t_under_cutoff=t_under_cutoff);\n  return(peaks);\n}\n\n\n\ndiff_peaks.plot <- function(peaks)\n{\n  x_label <- colnames(peaks$D)[peaks$ref_cols[1]];\n  y_label <- colnames(peaks$D)[peaks$signal_cols[1]];\n  x <- apply(peaks$D[,ref_cols],1,mean);\n  y <- apply(peaks$D[,signal_cols],1,mean);\n  xlim <- log2(c(min(c(x,y)),max(c(x,y))));\n  ylim <- log2(c(min(c(peaks$t_over,peaks$t_under)),max(c(peaks$t_over,peaks$t_under))));\n  smoothScatter(log2(x),log2(peaks$t_over),xlim=xlim,ylim=ylim,xlab=paste(x_label,' mean (log2)',sep=''),ylab='fold-change (log2)',main='signal vs reference fold-changes');\n  lines(log2(peaks$t_over_cutoff$value_bin),log2(peaks$t_over_cutoff$t_cutoff),col='red');\n  smoothScatter(log2(y),log2(peaks$t_under),xlim=xlim,ylim=ylim,xlab=paste(y_label,' mean (log2)',sep=''),ylab='fold-change (log2)',main='reference vs signal fold-changes');\n  lines(log2(peaks$t_under_cutoff$value_bin),log2(peaks$t_under_cutoff$t_cutoff),col='green');\n  ylim <- xlim;\n  smoothScatter(log2(x),log2(y),xlim=xlim,ylim=ylim,xlab=paste(x_label,' mean (log2)',sep=''),ylab=paste(y_label,' mean (log2)',sep=''),main='differential peaks');\n  pos_i <- peaks$t_over_cutoff$t_score>=1;\n  points(log2(x[pos_i]),log2(y[pos_i]),pch=19,col='red');\n  neg_i <- peaks$t_under_cutoff$t_score>=1;\n  points(log2(x[neg_i]),log2(y[neg_i]),pch=19,col='green');\n}\n\n\n\ndiff_peaks.store <- function(D,peaks,data_file)\n{\n  out_file <- paste(data_file,'.FDR=',fdr,'.score',sep='');\n  cat('output file = '); cat(out_file,fill=T); \n  out <- round(cbind(peaks$t_over_cutoff$t_score,peaks$t_under_cutoff$t_score,D),digits=6);\n  out <- cbind(rownames(D),out);\n  colnames(out)[1] <- 'locus';\n  colnames(out)[2] <- 'adjusted fold change (signal vs reference)';\n  colnames(out)[3] <- 'adjusted fold change (reference vs signal)';\n  colnames(out)[4] <- paste(colnames(out)[4],'normalized density 1')\n  colnames(out)[5] <- paste(colnames(out)[5],'normalized density 2')\n  colnames(out)[6] <- paste(colnames(out)[6],'normalized density 1')\n  colnames(out)[7] <- paste(colnames(out)[7],'normalized density 2')\n  write.table(out,out_file,quote=F,row.names=F,col.names=T,sep='\t');\n}\n\n\n\n\n\nremove_outliers <- function(D,outlier_pval)\n{\n  signal_z <- log2(D);\n  signal_label <- colnames(D)[1];\n  smoothScatter(signal_z,xlab=paste(signal_label,' replicate #1 (log2)',sep=''),ylab=paste(signal_label,' replicate #2 (log2)',sep=''),main='replicate reproducibility');\n  i <- seq(1,nrow(signal_z),length.out=20000);\n  fit <- loess(signal_z[i,2] ~ signal_z[i,1],span=0.5,degree=1);\n  x <- sort(signal_z[i,1]);\n  lines(x,predict(fit,x),col='magenta');\n  r <- signal_z[,2]-predict(fit,signal_z[,1]);\n  w <- dnorm(r,mean(r,na.rm=T),sd(r,na.rm=T))>outlier_pval;\n  w[is.na(w)] <- TRUE;\n  points(signal_z[!w,],pch=19,col='brown');\n  return(w);\n}\n\n\n\n\n\n\n#\n# MAIN PROGRAM\n#\n\n\nlibrary('preprocessCore');\nlibrary('MASS');\n\nargs <- commandArgs(trailingOnly=T);\ndata_file <- args[1];\nparam_file <- args[2];\nimage_file <- args[3];\nimage_type <- rev(unlist(strsplit(image_file,'.',fixed=T)))[1];\n\nparams <- readLines(param_file);\nn_signal_cols <- as.numeric(params[1]);\nn_ref_cols <- as.numeric(params[2]);\nwin_size <- as.numeric(params[3]);\nlogpval <- log(as.numeric(params[4]));\npseudo <- as.numeric(params[5]);		# add minimum possible pseudocount to avoid division by zero in fold-change computation\noutlier_pval <- as.numeric(params[6]);\nfdr <- as.numeric(params[7]);\nn_fdr_bins <- as.numeric(params[8]);\nsample_labels <- strsplit(params[9],',')[[1]];\nimage_size <- as.numeric(strsplit(params[10],',')[[1]]);\nimage_resolution <- as.numeric(params[11]);\n\n\n\ncat('Reading input matrix...',fill=T);\nDATA <- as.matrix(read.table(data_file,row.names=1,header=T,sep='\t',comment.char=''))[,1:(n_signal_cols+n_ref_cols)];\nif (nrow(DATA)==0) { cat('Error: data file is empty!',fill=T); q(); }\nif (n_signal_cols==1) {\n  DATA <- cbind(DATA[,1],sample_poisson(DATA[,1]),DATA[,2:ncol(DATA)]);\n  n_signal_cols <- 2;\n}\nif (n_ref_cols==1) {\n  DATA <- cbind(DATA[,1:(n_signal_cols+1)],sample_poisson(DATA[,n_signal_cols+1]));\n  n_ref_cols <- 2;\n}\nsignal_cols <- 1:n_signal_cols;\nref_cols <- n_signal_cols+1:n_ref_cols;\ncolnames(DATA) <- c(rep(sample_labels[1],n_signal_cols),rep(sample_labels[2],n_ref_cols));\n\n\n# compute diff peaks\nif (image_type=='tif') {\n  tiff(image_file,width=image_size[1],height=image_size[2],res=image_resolution,compression='lzw');\n} else if (image_type=='pdf') {\n  pdf(image_file); \n}\nlayout(matrix(1:6,ncol=3,byrow=F));\nD <- normalize_matrix((DATA+pseudo)/(win_size+pseudo));\nw_signal <- remove_outliers(D[,signal_cols],outlier_pval);\nw_ref <- remove_outliers(D[,ref_cols],outlier_pval);\npeaks <- diff_peaks.calc(D[w_signal&w_ref,],signal_cols,ref_cols,fdr,n_fdr_bins,win_size,logpval); \ndiff_peaks.plot(peaks);\ndiff_peaks.store(D[w_signal&w_ref,],peaks,data_file);\ndev.off();\n\n\n\ncat('Done.',fill=T);\ncat('*********************************************************************',fill=T);\n\n\n";

const char *RSCRIPT_TEMPLATE_PROFILE = \
"\n##\n## USAGE: genomic_apps.profile.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\nargs <- commandArgs(trailingOnly=T);\nprofile_data_file <- args[1];\nprofile_param_file <- args[2];\nprofile_image_file <- args[3];\nimage_type <- rev(unlist(strsplit(profile_image_file,'.',fixed=T)))[1];\n\nparams <- readLines(profile_param_file);\nprofile_bin_start <- as.numeric(params[1]);\nprofile_bin_stop <- as.numeric(params[2]);\nprofile_legend <- strsplit(params[3],',')[[1]];\nprofile_colors <- strsplit(params[4],',')[[1]];\nprofile_title <- params[5];\nprofile_xlab <- params[6];\nprofile_ylab <- params[7];\nimage_size <- as.numeric(strsplit(params[8],',')[[1]]);\nimage_resolution <- as.numeric(params[9]);\n  \nif (image_type=='tif') {\n  tiff(profile_image_file,width=image_size[1],height=image_size[2],res=image_resolution,compression='lzw');\n} else if (image_type=='pdf') {\n  pdf(profile_image_file); \n}\nY <- as.matrix(read.table(profile_data_file,header=FALSE,row.names=1,sep='\t'));\n\nd <- (profile_bin_stop-profile_bin_start)/ncol(Y);\nx <- seq(profile_bin_start+d/2,profile_bin_stop-d/2,d);\nplot(x,Y[1,],type='l',col=profile_colors[1],main=profile_title,xlab=profile_xlab,ylab=profile_ylab,xlim=c(profile_bin_start,profile_bin_stop),ylim=c(0,max(Y)));\ni <- 2;\nwhile (i <= nrow(Y)) \n{\n  lines(x,Y[i,],col=profile_colors[i]);\n  i <- i + 1;\n}\nlegend('topright',profile_legend,pch=16,col=profile_colors,inset=0.05);\ndev.off();\n\n  \n  \n\n\n";

