
library(tidyverse)
library(viridis)
library(plotly)
library(ggplot2)

#path to spreadsheet
NGS_cts_df <- data.frame(read.csv("/Users/nickles/Desktop/ALM_src_array_data.csv", header = TRUE))

#determine which rows to eliminate based on min initial reads
init_read_thresh <- 30
thresh_rows <- which(NGS_cts_df$X20_0 < init_read_thresh | NGS_cts_df$X20_10 < init_read_thresh | NGS_cts_df$X20_50 < init_read_thresh)
df_filt <- NGS_cts_df[-thresh_rows,]

num_plasmids <- length(df_filt$X)


#get means first
c12_00iptg_T0_tot <- sum(df_filt$X20_0)
c12_10iptg_T0_tot <- sum(df_filt$X20_10)
c12_50iptg_T0_tot <- sum(df_filt$X20_50)

c12_00iptg_G_tot <- sum(df_filt$X20_0.Top.5)
c12_10iptg_G_tot <- sum(df_filt$X20_10.Top.5)
c12_50iptg_G_tot <- sum(df_filt$X20_50.Top.5)


#now calc percent abundance at diff times and concentrations
pabund_00iptg_T0 <- df_filt$X20_0/rep(c12_00iptg_T0_tot, num_plasmids)
pabund_10iptg_T0 <- df_filt$X20_10/rep(c12_10iptg_T0_tot, num_plasmids)
pabund_50iptg_T0 <- df_filt$X20_50/rep(c12_50iptg_T0_tot, num_plasmids)

pabund_00iptg_growth <- df_filt$X20_0.Top.5/rep(c12_00iptg_G_tot, num_plasmids)
pabund_10iptg_growth <- df_filt$X20_10.Top.5/rep(c12_10iptg_G_tot, num_plasmids)
pabund_50iptg_growth <- df_filt$X20_50.Top.5./rep(c12_50iptg_G_tot, num_plasmids)

#enrichment
enrich_00iptg <- pabund_00iptg_growth/pabund_00iptg_T0
enrich_10iptg <- pabund_10iptg_growth/pabund_10iptg_T0
enrich_50iptg <- pabund_50iptg_growth/pabund_50iptg_T0

#make new df
df_enrich <- data.frame(cbind(df_filt$X, enrich_00iptg, enrich_10iptg, enrich_50iptg))


#eliminate more rows based on small change in enrichment bw iptg concentrations
enrich_fold_thresh <- 15
en_thresh_rows <- which(df_enrich$enrich_00iptg < enrich_fold_thresh & df_enrich$enrich_10iptg < enrich_fold_thresh & df_enrich$enrich_50iptg < enrich_fold_thresh)
df_enrich <- df_enrich[-en_thresh_rows,]
num_plasmids2 <- length(df_enrich$V1)

#make df for plotting
df2 <- data.frame(x = c(rep(0, num_plasmids2),rep(10,num_plasmids2), rep(50, num_plasmids2)), y = as.numeric(c(df_enrich$enrich_00iptg, df_enrich$enrich_10iptg, df_enrich$enrich_50iptg)), 
                  TS = as.character(rep(df_enrich$V1, 3)))


# ggp <- ggplot(df2, aes(x,y, col = TS)) +
#   geom_line(aes(color=TS), linewidth=1.5) +
#   geom_point(aes(color=TS), size=3.5) +
#   theme(legend.position="none") +
#   scale_colour_viridis(discrete = TRUE, option = "A", direction = -1, name = "[3-AT] (mM)" ) +
#   scale_fill_viridis(discrete = TRUE) +
#   ylim(0,15)
# 
# ggp







#ds <- t(data.frame(rbind(df_filt$X20_0.Top.5, df_filt$X20_10.Top.5)))

#df2 <- data.frame(x = c(rep(0, length(df_filt$X)),rep(10,length(df_filt$X)), rep(50, length(df_filt$X))), y = c(t(df_filt$X20_0.Top.5), t(df_filt$X20_10.Top.5), t(df_filt$X20_50.Top.5)), 
                 # TS = as.character(rep(df_filt$X, 3)))

#df2 <- data.frame(x = c(rep(0, length(df_filt$X)),rep(10,length(df_filt$X)), rep(50, length(df_filt$X))),
                #  TS = as.character(rep(df_filt$X, 3)))


#plot using plotly


#df2 <- rbind(enrich_00iptg, enrich_10iptg, enrich_50iptg)
#x <- c(rep(c(0, 10, 50), num_plasmids2))

#ptly_df <- data.frame(x, df2)

p <- plot_ly() %>% layout(showlegend = FALSE)

for(i in 1:num_plasmids2){
    p <- add_trace(p, x = c(0,10,50), y = as.numeric(df_enrich[i,2:4]), type='scatter', mode='line', name = df_enrich$V1[i])
}

p

########################

# 
# 
# ggp <- ggplot(df2, aes(x,y, col = TS)) +
#   geom_line(aes(color=TS), linewidth=1.5) +
#   geom_point(aes(color=TS), size=3.5) +
#   theme(legend.position="none") +
#   scale_colour_viridis(discrete = TRUE, option = "A", direction = -1, name = "[3-AT] (mM)" ) +
#   scale_fill_viridis(discrete = TRUE) +
#   ylim(0,1000)
# 
# ggp
# 
# 
# 
# 
# 
# 
# #df_filt <- data.frame(read.csv("/Users/nickles/Desktop/NGS/NGS_cts_test1.csv", header = TRUE, nrows = 5))
# 
# df_filt <- setNames(data.frame(t(df_filt[,-1])), paste0(df_filt[,1]))
# 
# rownames(df_filt) <- substr(rownames(df_filt), 1, nchar(rownames(df_filt))-1) # remove . in names
# 
# Areps <- c(paste0(A_cond, '_R1'), paste0(A_cond, '_R2'))
# Breps <- c(paste0(B_cond, '_R1'), paste0(B_cond, '_R2'))
# 
# # removal of members prior to deseq?
# NGS_cts_filt <- subset(df_filt, Areps[1] > read_thresh & Areps[2] > read_thresh & 
#                          Breps[1] > read_thresh & Breps[2] > read_thresh) #filter plasmids with few reads
# 
# countdata <- as.matrix(bind_cols(NGS_cts_filt[Areps[1]], NGS_cts_filt[Areps[2]], NGS_cts_filt[Breps[1]], NGS_cts_filt[Breps[2]]))
# countdata <- as.matrix(bind_cols(df_filt[Areps[1]], NGS_cts_filt[Areps[2]], NGS_cts_filt[Breps[1]], NGS_cts_filt[Breps[2]]))
