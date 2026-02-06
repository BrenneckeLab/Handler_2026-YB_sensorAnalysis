usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    devtools::install_github("tidyverse/ggplot2")
  require(p, character.only = TRUE)
}
usePackage("ggplot2")

usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, ,repos = "http://cran.us.r-project.org")
  require(p, character.only = TRUE)
}

usePackage("dplyr")
usePackage("tibble")
usePackage("tidyr")
usePackage("readr")
usePackage("stringr")
usePackage("cowplot")
usePackage("plotly")
usePackage("RColorBrewer")

theme_set(theme_cowplot())

###################################################################################################
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################
#TMP="~/Desktop/files for eln/"
setwd( TMP)

#NAME="empty_biogenSensor_total"

fullNAME=paste(TMP,NAME,"_sense.norm.bg", sep="")
RAWplus=read_tsv(fullNAME)
# fullNAME=paste(TMP,NAME,"_antisense.bg", sep="")
# RAWminus=read_tsv(fullNAME, col_names = c("sensor", "coord", "value_minus"), cols(value_minus=col_double()))
safeNAME=gsub("-",".",NAME)

#determine if both strands are analyzed or not
#RAW = bind_cols(RAWplus, RAWminus)
RAW=RAWplus

#read in coordinates for sensor parts
fullNAME=paste(TMP,"sensor_parts.bed", sep="")
ANN = read_tsv(fullNAME, col_names = c("SENSOR", "START", "END", "NAME", "X", "STRAND"))
ANN

# #prepare normalization factors
# #sequencing depth normalization
# miRNA_NORM=as.numeric(miRNA_NORM)
# piRNA_NORM=as.numeric(piRNA_NORM)

# #sensor normalization
# shRNA_NORM=as.numeric(shRNA_NORM)
# shRNA_NORM=shRNA_NORM/miRNA_NORM #shRNA norm needs to be corrected for sequencing depth
# QPCR_NORM=as.numeric(QPCR_NORM)

# #normalize to miRNA, shRNA and piRNA 
# RAW = RAW %>%
#   mutate(
#     sRNA_miRNA = sRNA_raw/miRNA_NORM,
#     sRNA_shRNA = sRNA_miRNA/shRNA_NORM,
#     sRNA_QPCR = sRNA_shRNA/QPCR_NORM
#   )


# maxVAL= RAW %>%
#   filter( coord>3000, sensor != "tj-locus" ) %>%
#   summarise(xy=max(value_plus_shmiR))%>%
#   pull(xy)

RAW 

ANN = ANN %>%
  separate(SENSOR, into=c("Plasmid","SENSOR"), sep="_", extra="merge")
ANN

RAW
#shRNA normalized
p = RAW %>% 
    pivot_longer(cols=c(COUNT_shNORM, COUNT_QPCRnorm, COUNT_miRNAnorm), names_to="NORMtype", values_to="VALUE")%>%
    ggplot() + 
      # geom_rect(data=ANN, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), alpha=0.2, colour="gray60", fill="gray80") +
      geom_line(aes(x=POSITION, y=VALUE)) +
      # geom_text(data=ANN, aes(x=START+(END-START)/2, y=-5, label=NAME), size=2) +
      facet_grid(SENSOR~NORMtype, scales="free_y", labeller =label_wrap_gen(width = 50, multi_line = TRUE) )+
      labs(title=NAME)+
      theme(strip.text = element_text(size = 5))
      {}

FILENAME=paste(OPENdir, NAME, ".png", sep="")
ggsave(FILENAME,p, width=30, height=15, dpi=300 ,units="cm")

# #QPCR normalized
# p = RAW %>% 
#       # geom_rect(data=ANN, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), alpha=0.2, colour="gray60", fill="gray80") +
#       geom_line(aes(x=POSITION, y=sRNA_QPCR)) +
#       # facet_wrap(~SENSOR, scales="free_y", ncol=1, labeller =label_wrap_gen(width = 50, multi_line = TRUE) )+
#       labs(title=NAME)+
#       theme(strip.text = element_text(size = 5))
#       {}

# FILENAME=paste(OPENdir, NAME, ".QPCR.png", sep="")
# ggsave(FILENAME,p, width=30, height=15, dpi=300 ,units="cm")

q()
EXPORT= filter(RAW, sensor == "tj-locus") %>%
  select(value_plus)
# EXPORT=as.dataframe(EXPORT)
EXPORT=rename(EXPORT, !!NAME:=value_plus)
EXPORT
FILENAME=paste(RAWdir, NAME,".tj-locus.txt", sep="")
write_delim(EXPORT,FILENAME)


EXPORT= filter(RAW, !grepl("tj-locus",sensor)) %>%
  filter(!grepl("myc",sensor)) %>%
  mutate(ID = NAME) %>%
  select(sensor, coord,value_plus_shmiR,value_plus, ID)


# EXPORT=as.dataframe(EXPORT)
#EXPORT=rename(EXPORT, !!NAME:=value_plus)
EXPORT
FILENAME=paste(RAWdir, NAME,".txt", sep="")
write_delim(EXPORT,FILENAME, col_names=FALSE)

# write.table(EXPORT,file=FILENAME, append=TRUE)

# pp = ggplotly(p)
# # FILENAME=paste(OPENdir, NAME, ".html", sep="")
# #
#  htmlwidgets::saveWidget(pp, FILENAME)
# RAWzoom=filter(RAW,coord>=3575 & coord<=9000)
# ANNzoom=filter(ANN,start>=3575 & start<=9000)
# RAWzoom

# p = ggplot() + 
#       geom_line(data=filter(RAWzoom, grepl('sBR',sensor)), aes(x=coord, y=value_plus)) +
#       geom_line(data=filter(RAWzoom, grepl('sBR',sensor)), aes(x=coord, y=value_minus)) +
#       #scale_y_continuous(limits=c(-100,1200))+
#       scale_x_continuous(limits=c(3575,9000))+
#       geom_rect(data=filter(ANNzoom, grepl('sBR',sensor)), aes(xmin=start, xmax=end, ymin=-25, ymax=-50),  colour = "gray60", fill = "gray80")+
#       geom_text(data=filter(ANNzoom, grepl('sBR',sensor)), aes(x=start+(end-start)/2, y=-90, label=name, colour=), size=2) +
#       {}

# FILENAME=paste(OPENdir, NAME, ".zoom.pdf", sep="")
# ggsave(FILENAME,p, width=30, height=15, dpi=300 ,units="cm")
# RAW <- rename(RAW, !!safeNAME := "value_plus")
# head(RAW)
# shmirNAME=paste(safeNAME,".shmiR",sep="")
# RAW <- rename(RAW, !!shmirNAME := "value_plus_shmiR")

# head(RAW)
# colnames(RAW)[colnames(RAW) == "value_plus_shmiR"]<paste(NAME,".shmiR",sep="")
# FILENAME=paste(OPENdir, NAME, ".table.PLH.txt", sep="")
# write.table(filter(RAW, grepl('PLH',sensor)),FILENAME,row.names=FALSE,sep="\t", quote = FALSE)
# FILENAME=paste(OPENdir, NAME, ".table.tj.txt", sep="")
# write.table(filter(RAW, grepl('tj',sensor)),FILENAME,row.names=FALSE,sep="\t", quote = FALSE)

# shCOUNT=RAW[2591,"value_plus"]
# tjARRAY=filter(RAW, grepl('tj',sensor))
# tjVALUES=tjARRAY[3100:4500,"value_plus"]
# tjCOUNT=median(tjVALUES[[1]])
# RATIO=shCOUNT/tjCOUNT

# x=c(NAME,NORM,tjCOUNT,shCOUNT,RATIO)
# y=unlist(x, use.names = FALSE,)
# FILENAME=paste(OPENdir, "stats.tmp", sep="")
# write(paste(y,collapse=" "),file=FILENAME, append=TRUE)