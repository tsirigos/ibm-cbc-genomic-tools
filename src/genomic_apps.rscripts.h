const char *RSCRIPT_TEMPLATE_HEATMAP = \
"\n##\n## USAGE: genomic_apps.heatmap.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\nnorm_rows <- function(X) { for (i in 1:nrow(X)) X[i,] <- (X[i,]-mean(X[i,]))/sd(X[i,]); return(X); }\n\nargs <- commandArgs(trailingOnly=T);\ndata_file <- args[1];\nparam_file <- args[2];\nimage_file <- args[3];\n\nparams <- readLines(param_file);\nshift_upstream <- as.numeric(params[1]);\nshift_downstream <- as.numeric(params[2]);\nheatmap_colors <- strsplit(params[3],',')[[1]];\nheatmap_title <- strsplit(params[4],',')[[1]];\nheatmap_xlab <- params[5];\nheatmap_ylab <- params[6];\nheatmap_size <- as.numeric(strsplit(params[7],',')[[1]]);\nheatmap_resolution <- as.numeric(params[8]);\nn_heatmaps <- as.numeric(params[9]);\n\n\nlibrary('MASS');\nlibrary('preprocessCore');           # from Bioconductor\nlibrary('gplots');                   # from Bioconductor\n\n# load data\nD <- as.matrix(read.table(data_file,row.names=1,sep='\t'));\n  \n# create combined heatmap (main version)\ntiff(image_file,width=heatmap_size[1],height=heatmap_size[2],res=heatmap_resolution,compression='lzw');\n\npar(fig=c(0,1,0,1),mar=c(2,2,0,0)); \nplot.new();\nmtext(heatmap_xlab,side=1);\nmtext(heatmap_ylab,side=2);\n\n\nd <- ncol(D)/n_heatmaps;\nI <- 1:d;\ndj <- 0.9/n_heatmaps;\nfor (j in 1:n_heatmaps) {\n  par(fig=c(0.1+(j-1)*dj,0.1+j*dj,0.1,1),mar=c(0.5,0.5,3,0.5),new=TRUE);\n  colorscale <- c(colorpanel(20,low='white',high='white'),colorpanel(50,low='white',high=heatmap_colors[j]),colorpanel(30,low=heatmap_colors[j],high=heatmap_colors[j]));\n  image(z=t(norm_rows(D[,I])),col=colorscale,main=heatmap_title[j],xlab=heatmap_xlab,ylab=heatmap_ylab,xaxt='n',yaxt='n');\n  I <- I+d;\n}\n\ndev.off();\n\n";

const char *RSCRIPT_TEMPLATE_PEAKDIFF = \
"# source('c:/Aris/Research/Code/genomic_tools/genomic_apps.peakdiff.r');\n\n\n##\n## USAGE: genomic_apps.peakdiff.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\n\nmax_abs <- function(a,b) \n{\n  if (abs(a)>abs(b)) a else b\n}\n\nscore_binomial <- function(x,y,n) \n{\n  a <- -pbinom(n*x,n,y,lower.tail=F,log.p=T);\n  b <- pbinom(n*x,n,y,log.p=T);\n  max_abs(a,b)\n}\n\nscore_fold <- function(x,y) \n{\n  return(x/y);\n}\n\nscore <- function(x,y,n,logpval) \n{\n  p <- score_binomial(x,y,n);\n  return(if (abs(p) < -logpval) 1.0 else score_fold(x,y));\n}\n\n\ncalc_fdr_cutoff <- function(pos,neg,fdr) \n{\n  if (fdr<=0) return(Inf);\n  kpos <- 1;\n  kneg <- 1;\n  while ((kpos<=length(pos))&(kneg<=length(neg))) {\n    if ((length(neg)-kneg+1)/(length(pos)-kpos+1)<=fdr) { break; }\n    if (pos[kpos]<neg[kneg]) { kpos <- kpos+1; }\n    else if (pos[kpos]>neg[kneg]) { kneg <- kneg+1; }\n    else { kpos <- kpos+1; kneg <- kneg+1; }\n  }\n  if (kpos>length(pos)) { y <- Inf; }\n  else { y <- pos[kpos]; }\n  return(y);\n}\n\n\nmean_paired_score <- function(signal,ref,win_size,logpval)\n{\n  for (i in 1:length(signal)) s <- s + score(signal[i],ref[i],win_size,logpval);\n  return(s/length(signal));\n}\n\n\npaired_scores <- function(signal,ref,win_size,logpval)\n{\n  s <- signal;\n  for (i in 1:length(signal)) s[i] <- score(signal[i],ref[i],win_size,logpval);\n  return(s);\n}\n\n\ncalc_diff_peaks.plain <- function(D,signal_cols,ref_cols,fdr,win_size,logpval)\n{\n  cat('Computing scores...\n');\n  t_over <- apply(D,1,function(x) score(mean(x[signal_cols]),mean(x[ref_cols]),win_size,logpval));\n  t_under <- apply(D,1,function(x) score(mean(x[ref_cols]),mean(x[signal_cols]),win_size,logpval));\n  #t_over <- apply(D,1,function(x) min(paired_scores(x[signal_cols],x[ref_cols],win_size,logpval)));      # use these for paired samples\n  #t_under <- apply(D,1,function(x) min(paired_scores(x[ref_cols],x[signal_cols],win_size,logpval)));\n  cat('Computing randomized scores...\n');\n  t_control <- apply(D,1,function(x) mean(c(score(x[signal_cols[1]],x[signal_cols[2]],win_size,logpval),score(x[ref_cols[1]],x[ref_cols[2]],win_size,logpval),score(x[signal_cols[2]],x[signal_cols[1]],win_size,logpval),score(x[ref_cols[2]],x[ref_cols[1]],win_size,logpval))));\n  cat('Computing FDR scores...\n');\n  t_over_cutoff <- as.vector(calc_fdr_cutoff(sort(t_over),sort(t_control),fdr));\n  t_under_cutoff <- as.vector(calc_fdr_cutoff(sort(t_under),sort(t_control),fdr));\n  t_over_score <- cbind(t_over,D)[(t_over>=t_over_cutoff),];\n  t_under_score <- -cbind(t_under,D)[(t_under>=t_under_cutoff),];\n  cat('significant entries (over) = '); cat(nrow(t_over_score)); cat('\n');\n  cat('significant entries (under) = '); cat(nrow(t_under_score)); cat('\n');\n  out <- list(D=D,signal_cols=signal_cols,ref_cols=ref_cols,fdr=fdr,t_over=t_over,t_under=t_under,t_control=t_control,t_over_cutoff=t_over_cutoff,t_under_cutoff=t_under_cutoff,t_over_score=t_over_score,t_under_score=t_under_score);\n  calc_diff_peaks.plain.plot.mean(out);\n  return(out);\n}\n\n\n\ncalc_diff_peaks.plain.plot.mean <- function(out)\n{\n  x_label <- colnames(out$D)[out$ref_cols[1]];\n  y_label <- colnames(out$D)[out$signal_cols[1]];\n  x <- apply(out$D[,ref_cols],1,mean);\n  y <- apply(out$D[,signal_cols],1,mean);\n  xlim <- ylim <- log(c(min(c(x,y)),max(c(x,y))));\n  smoothScatter(log(x),log(y),xlim=xlim,ylim=ylim,xlab=paste(x_label,' mean (log)',sep=''),ylab=paste(y_label,' mean (log)',sep=''),main='step 2: differential peak detection');\n  pos_i <- out$t_over>=out$t_over_cutoff;\n  points(log(x[pos_i]),log(y[pos_i]),pch=19,col='red');\n  neg_i <- out$t_under>=out$t_under_cutoff;\n  points(log(x[neg_i]),log(y[neg_i]),pch=19,col='green');\n}\n\n\n\n\ncalc_diff_peaks.plain.plot <- function(out)\n{\n  x_label <- colnames(out$D)[out$ref_cols[1]];\n  y_label <- colnames(out$D)[out$signal_cols[1]];\n  x <- out$D[,ref_cols[1]];\n  y <- out$D[,signal_cols[1]];\n  smoothScatter(log(x),log(y),xlab=paste(x_label,' replicate #1 (log)',sep=''),ylab=paste(y_label,' replicate #1 (log)',sep=''),main='step 2: differential peak detection');\n  pos_i <- out$t_over>=out$t_over_cutoff;\n  points(log(x[pos_i]),log(y[pos_i]),pch=19,col='red');\n  neg_i <- out$t_under>=out$t_under_cutoff;\n  points(log(x[neg_i]),log(y[neg_i]),pch=19,col='green');\n  x <- out$D[,ref_cols[2]];\n  y <- out$D[,signal_cols[2]];\n  smoothScatter(log(x),log(y),xlab=paste(x_label,' replicate #2 (log)',sep=''),ylab=paste(y_label,' replicate #2 (log)',sep=''),main='step 2: differential peak detection');\n  pos_i <- out$t_over>=out$t_over_cutoff;\n  points(log(x[pos_i]),log(y[pos_i]),pch=19,col='red');\n  neg_i <- out$t_under>=out$t_under_cutoff;\n  points(log(x[neg_i]),log(y[neg_i]),pch=19,col='green');\n}\n\n\n\n\n\nremove_outliers <- function(D,outlier_pval)\n{\n  signal_z <- log(D);\n  signal_label <- colnames(D)[1];\n  smoothScatter(signal_z,xlab=paste(signal_label,' replicate #1 (log)',sep=''),ylab=paste(signal_label,' replicate #2 (log)',sep=''),main='step 1: replicate reproducibility');\n  i <- seq(1,nrow(signal_z),length.out=20000);\n  fit <- loess(signal_z[i,2] ~ signal_z[i,1],span=0.5,degree=1);\n  x <- sort(signal_z[i,1]);\n  lines(x,predict(fit,x),col='magenta');\n  r <- signal_z[,2]-predict(fit,signal_z[,1]);\n  w <- dnorm(r,mean(r,na.rm=T),sd(r,na.rm=T))>outlier_pval;\n  w[is.na(w)] <- TRUE;\n  points(signal_z[!w,],pch=19,col='brown');\n  return(w);\n}\n\n\n\nnormalize_matrix <- function(D)\n{\n  D_norm <- normalize.quantiles(D);\n  dimnames(D_norm) <- dimnames(D);\n  return(D_norm);\n}\n\n\n\n\n\n\n#\n# MAIN PROGRAM\n#\n\n\nlibrary('preprocessCore');\nlibrary('MASS');\n\nargs <- commandArgs(trailingOnly=T);\ndata_file <- args[1];\nparam_file <- args[2];\nimage_file <- args[3];\n\nparams <- readLines(param_file);\nn_signal_cols <- as.numeric(params[1]);\nn_ref_cols <- as.numeric(params[2]);\nwin_size <- as.numeric(params[3]);\nlogpval <- log(as.numeric(params[4]));\npseudo <- as.numeric(params[5]);		# add minimum possible pseudocount to avoid division by zero in fold-change computation\noutlier_pval <- as.numeric(params[6]);\nfdr <- as.numeric(params[7]);\nsample_labels <- strsplit(params[8],',')[[1]];\nimage_size <- as.numeric(strsplit(params[9],',')[[1]]);\nimage_resolution <- as.numeric(params[10]);\n\n\n\n\ncat('Reading input matrix...\n');\nDATA <- as.matrix(read.table(data_file,row.names=1,sep='\t'));\nif (nrow(DATA)==0) { cat('Error: data file is empty!\n'); q(); }\n#DATA <- DATA[seq(1,nrow(DATA),length.out=100000),];\nsignal_cols <- 1:n_signal_cols;\nref_cols <- n_signal_cols+1:n_ref_cols;\ndimnames(DATA)[[2]][signal_cols] = sample_labels[1];\ndimnames(DATA)[[2]][ref_cols] = sample_labels[2];\n\n\n# compute diff peaks\ntiff(image_file,width=image_size[1],height=image_size[2],res=image_resolution,compression='lzw');\nlayout(matrix(1:4,ncol=2,byrow=TRUE));\nD <- normalize_matrix((DATA+pseudo)/(win_size+pseudo));\nw_signal <- remove_outliers(D[,signal_cols],outlier_pval);\nw_ref <- remove_outliers(D[,ref_cols],outlier_pval);\nout <- calc_diff_peaks.plain(D[w_signal&w_ref,],signal_cols,ref_cols,fdr,win_size,logpval); \ndev.off();\nout_file <- paste(data_file,'.FDR=',fdr,'.score',sep='');\ncat('output file = '); cat(out_file); cat('\n');\nwrite.table(format(rbind(out$t_over_score,out$t_under_score),digits=4),out_file,quote=F,col.names=F,sep='\t');\n\n\n\ncat('Done.\n');\ncat('*********************************************************************\n');\n\n\n";

const char *RSCRIPT_TEMPLATE_PROFILE = \
"\n##\n## USAGE: genomic_apps.profile.r DATA-FILE PARAMETER-FILE OUTPUT-IMAGE-FILE\n##\n\nargs <- commandArgs(trailingOnly=T);\nprofile_data_file <- args[1];\nprofile_param_file <- args[2];\nprofile_image_file <- args[3];\n\nparams <- readLines(profile_param_file);\nprofile_upstream <- as.numeric(params[1]);\nprofile_downstream <- as.numeric(params[2]);\nprofile_legend <- strsplit(params[3],',')[[1]];\nprofile_colors <- strsplit(params[4],',')[[1]];\nprofile_title <- params[5];\nprofile_xlab <- params[6];\nprofile_ylab <- params[7];\n  \ntiff(profile_image_file,width=3500,height=3500,res=600,compression='lzw',antialias='none');\nY <- as.matrix(read.table(profile_data_file,header=FALSE,row.names=1,sep='\t'));\n\nd <- (profile_upstream+profile_downstream)/ncol(Y);\nx <- seq(-profile_upstream+d/2,+profile_downstream-d/2,d);\nplot(x,Y[1,],type='l',col=profile_colors[1],main=profile_title,xlab=profile_xlab,ylab=profile_ylab,xlim=c(-profile_upstream,profile_upstream),ylim=c(0,max(Y)));\ni <- 2;\nwhile (i <= nrow(Y)) \n{\n  lines(x,Y[i,],col=profile_colors[i]);\n  i <- i + 1;\n}\nlegend('topright',profile_legend,pch=16,col=profile_colors,inset=0.05);\ndev.off();\n\n  \n  \n\n\n";

