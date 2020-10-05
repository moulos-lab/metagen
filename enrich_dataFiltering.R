library(readr)
library(scales)

# Set script args
args = commandArgs(trailingOnly=TRUE)
path   <- args[1]

# Set config vars
fastq_path      <- paste0(path,"/fastq")
out_path        <- paste0(path,"/out/")
track           <- read.csv(paste0(out_path,"/data_trunc.tsv"), sep="")
sampleNames     <- rownames(track)

# for (i in sampleNames) {
# mapings <- read_delim(paste0("out/plots/reads_distributions/",i,"_filt_screen.txt"), 
#                       "\t", 
#                       escape_double = FALSE, 
#                       col_names = TRUE, 
#                       trim_ws = TRUE, 
#                       skip = 1
# )
# 
# mapings <- mapings[1:2,1:3]
# track[i,7] <- mapings[1,3]
# track[i,8] <- mapings[2,3]
# }

match_per <- function(sampleName) {
mapings <- read_delim(paste0(out_path,"/plots/reads_distributions/",sampleName,"_filt_screen.txt"),
                      "\t",
                      escape_double = FALSE,
                      col_names = TRUE,
                      trim_ws = TRUE,
                      skip = 1
)

mapings <- mapings[1:2,1:3]
track[sampleName,7] <- mapings[1,3]
track[sampleName,8] <- mapings[2,3]
return(track)  
}

track <- do.call(rbind,lapply(sampleNames, match_per))

track <- track[complete.cases(track),]
row.names(track) <- sampleNames
colnames(track) <- c("Total Input","QC Filtered","Denoised","Chimeras Filtered","#Bacteria", "%Bacteria", "#Human", "#Mouse")

track$`#Human` <- track$`QC Filtered` - track$`#Human`
track$`#Mouse` <- track$`QC Filtered` - track$`#Mouse`

track$`%Human` <- track$`#Human` / track$`QC Filtered`
track$`%Mouse` <-track$`#Mouse` / track$`QC Filtered`

track <- track[, c("Total Input","QC Filtered","Denoised","Chimeras Filtered","#Bacteria", "%Bacteria", "#Human", "%Human", "#Mouse", "%Mouse")]

track$`%Bacteria` <- label_percent(accuracy = 0.1)(track[,c("%Bacteria")])
track$`%Human`    <- label_percent(accuracy = 0.1)(track[,c("%Human")])
track$`%Mouse`    <- label_percent(accuracy = 0.1)(track[,c("%Mouse")])

write.table(track, file = paste0(out_path,"data_trunc.tsv"))
