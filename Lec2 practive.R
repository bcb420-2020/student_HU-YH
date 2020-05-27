#r, eval=FALSE}
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()


#r practive lecture code, echo=FALSE}
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
geo_tables <- dbListTables(con)
geo_tables
dbListFields(con,'gse')
results <- dbGetQuery(con,'select * from gpl limit 5')
#knitr::kable(head(results[,1:5]), format = "html")
num_platforms <- dbGetQuery(con,'select count(*) from gpl')
num_platforms
dbListFields(con,'gpl')
uniq_tech <- dbGetQuery(con,'select distinct technology from gpl')
uniq_tech
nrow(uniq_tech)
num_uniq_tech <- dbGetQuery(con,'select technology,count(*) from gpl group by technology')
colnames(num_uniq_tech)[2] <- "Num_Platforms"
plot_df <- num_uniq_tech[!is.na(num_uniq_tech$technology),]
p<-ggplot(data=plot_df, aes(technology, Num_Platforms)) +
  geom_col() + coord_flip()
p

species_ids <- dbGetQuery(con,'select organism,count(*) as num_plat from gpl where organism like "%homo%" group by organism order by num_plat desc')
knitr::kable(species_ids[1:5,], format = "html")

num_uniq_tech_human <- dbGetQuery(con,'select technology,count(*) as num_plat from gpl where organism = "Homo sapiens" group by technology  order by num_plat desc')
colnames(num_uniq_tech_human)[2] <- "Num_Platforms"

# r search for entry, eval = TRUE}
sql <- paste("SELECT DISTINCT gse.title,gse.gse, gpl.title,",
             " gse.submission_date",
             "FROM",
             "  gse JOIN gse_gpl ON gse_gpl.gse=gse.gse",
             "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
             "WHERE",
             "  gse.submission_date > '2016-01-01' AND",
             "  gse.title LIKE '%treatment%' AND",
             "  gpl.organism LIKE '%Homo sapiens%' AND",
             "  gpl.title LIKE '%HiSeq%' ",sep=" ")

rs <- dbGetQuery(con,sql)
dim(rs)

unlist(lapply(rs$supplementary_file,
              FUN = function(x){x <- unlist(strsplit(x,";")) ;
              x <- x[grep(x,pattern="txt",ignore.case = TRUE)];
              tail(unlist(strsplit(x,"/")),n=1)})) [1:10]
rs <- dbGetQuery(con,sql)
counts_files <- rs$supplementary_file[grep(rs$supplementary_file,
                                           pattern = "count",ignore.case = TRUE)]
dim(rs)
knitr::kable(rs, format = "html")
typeof(rs$submission_date)
sfiles = getGEOSuppFiles('GSE70072')
fnames = rownames(sfiles)
# there is only one supplemental file
b2 = read.delim(fnames[1],header=TRUE)
head(b2)

gse <- getGEO("GSE102004",GSEMatrix=FALSE)
kable(data.frame(head(Meta(gse))), format = "html")

C = read.delim('G:/BCB420Project/GSE70072/GSE70072_HumanSerousCancer_rawCounts.txt.gz', header=TRUE, check.names = FALSE)

85001 100860 104858
sample
106169
108539
109720
113493
considerate
110536



counts_density <- apply(log2(cpm(data_filtered[,2:7])), 2, density)
#calculate the limits across all the samples
xlim <- 0; ylim <- 0
for (i in 1:length(counts_density)) {
  xlim <- range(c(xlim, counts_density[[i]]$x)); 
  ylim <- range(c(ylim, counts_density[[i]]$y))
}
cols <- rainbow(length(counts_density))
ltys <- rep(1, length(counts_density))
#plot the first density plot to initialize the plot
plot(counts_density[[1]], xlim=xlim, ylim=ylim, type="n", 
     ylab="Smoothing density of log2-CPM", main="", cex.lab = 0.85)
#plot each line
for (i in 1:length(counts_density)) lines(counts_density[[i]], col=cols[i], lty=ltys[i])
#create legend
legend("topright", colnames(data2plot),  
       col=cols, lty=ltys, cex=0.75, 
       border ="blue",  text.col = "green4", 
       merge = TRUE, bg = "gray90")

plotMeanVar(d, show.raw.vars = TRUE,
            show.tagwise.vars=FALSE, NBline=FALSE, 
            show.ave.raw.vars = FALSE,show.binned.common.disp.vars = FALSE)


data2plot <- log2(data_filtered[,2:7])
boxplot(data_filtered[,-1], xlab = "Samples", ylab = "log2 CPM", 
        las = 2, cex = 0.5, cex.lab = 0.5,
        cex.axis = 0.5, main = "C4-12ERaERE RNASeq Samples")

dgList <- calcNormFactors(filtered_data_matrix, method="TMM")
