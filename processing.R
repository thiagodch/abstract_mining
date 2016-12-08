library("tm")
library("SnowballC")
library("ggfortify")
library("ggplot2")
library("wordcloud")
library("cluster")
library("proxy")
library("WGCNA")
library("corrplot")
library("igraph")
#library("topicmodels")

# Get abstracts and separate them into files
# esearch -db pubmed -query '"juvenile idiopathic arthritis"' | efetch -format medline > abstracts2.txt
# perl get_ab.pl > abs_table.txt

# Create corpus from abstract files
docs <- Corpus(DirSource("./abstracts_final"))

# Check contents of the 30th text file (example)
writeLines(as.character(docs[[30]]))

# Check lengths of abstracts (number of characters)
lens <- sapply(docs, function(x) nchar(x$content))
hist(unlist(lens), xlab="Tamanho", ylab="Frequência", main="")

# Pre-processing
## Transform colons, hyphens and other things into space
toSpace <- content_transformer(function(x, pattern) {return (gsub(pattern, " ", x))})
docs <- tm_map(docs, toSpace, "-")
docs <- tm_map(docs, toSpace, ":")
docs <- tm_map(docs, toSpace, "\\*")
docs <- tm_map(docs, toSpace, "'")
docs <- tm_map(docs, toSpace, " -")
docs <- tm_map(docs, toSpace, "\\.")
docs <- tm_map(docs, toSpace, '"')
## Replace proper punctuation signs with space
docs <- tm_map(docs, removePunctuation)

## Transform everything to lower case 
docs <- tm_map(docs, content_transformer(tolower))

## Remove numbers
docs <- tm_map(docs, removeNumbers)

## Remove stop words (a, an, the...)
docs <- tm_map(docs, removeWords, stopwords("english"))

## Remove extraneous whitespaces
docs <- tm_map(docs, stripWhitespace)

# Stemming
## Reduce words to their common stems (might not be so good?)
docs <- tm_map(docs, stemDocument)

## Remove additional stop words
myStopwords <- c("although", "use", "also", "can", "say", "one", "why", "use", "also", 
                 "howev", "tell", "will", "much", "need", "take", 
                 "tend", "even", "like", "particular", "rather", "said",
                 "get", "well", "make", "ask", "come", "end", 
                 "first", "two", "help", "often", "may", "might", 
                 "see", "someth", "thing", "point", "post", "look",
                 "right", "now", "think", "'ve", "'re", "s", "anoth","put","set","new","good",
                 "want","sure","kind","larg","yes,","day","etc",
                 "quit","sinc","attempt","lack","seen","awar",
                 "littl","ever","moreov","though","found","abl",
                 "enough","far","earli","away","achiev","draw",
                 "last","never","brief","bit","entir","brief",
                 "great","lot")
docs <- tm_map(docs, removeWords, myStopwords)

## Remove extraneous whitespaces
docs <- tm_map(docs, stripWhitespace)

## Fix bad stemmings
docs <- tm_map(docs, content_transformer(gsub), pattern="juvenileadult", replacement="juvenil")
docs <- tm_map(docs, content_transformer(gsub), pattern="juvennil", replacement="juvenil")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthralgi", replacement="arthralgia")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthralgiaarthr", replacement="arthralgia")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthralgiasarthr", replacement="arthralgia")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthrit", replacement="arthriti")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthrithi", replacement="arthriti")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthritid", replacement="arthriti")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthritidesspondyloarthropathi", replacement="arthriti")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthritidi", replacement="arthriti")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthritii", replacement="arthriti")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthritiogen", replacement="arthriti")
docs <- tm_map(docs, content_transformer(gsub), pattern="arthritisss", replacement="arthriti")

## Remove extraneous whitespaces
docs <- tm_map(docs, stripWhitespace)

# The document term matrix
## Matrix with all occurrences of words in the corpus, by document
## Rows: documents; Columns: words/terms
dtm <- DocumentTermMatrix(docs)
## Info
dtm
## Inspect dtm
inspect(dtm[1:2,1140:1150])

# Mining the corpus

## Getting the frequency of occurence of each word in the corpus
freq <- colSums(as.matrix(dtm))
## Create sort order (descending)
ord <- order(freq, decreasing = TRUE)
## Inspect most and least frequently occurring terms
freq[head(ord)]
freq[tail(ord)]

## Words with low frequencies tend to be more descriptive of specific documents
### Use inverse document frequencies to balance between frequency and specificity
### Remove words that give few information (like "can" and "one"). 
(dtmr <- DocumentTermMatrix(docs, control = list(weighting = function(x) weightTfIdf(x, normalize = FALSE), 
                                          wordLengths=c(4, 20), bounds=list(global = c(347, 3128)))))
### Only words that occur in at least 347 and at most 3128 documents 
### (10% and 90%) and which are between 4 and 20 characters in length

## Recalculate frequencies
freqr <- colSums(as.matrix(dtmr))
## Lengths should be total number of terms
length(freqr)
## Create sort order (desc)
ordr <- order(freqr, decreasing=TRUE)
## Inspect most and least frequently occurring terms
freqr[head(ordr)]
freqr[tail(ordr)]

## Find most frequent terms (at least 500 times in the corpus)
findFreqTerms(dtmr, lowfreq = 500)

## Find correlations between these and other terms in the corpus
### Quantitative measure of co-occurrence of words in multiple documents
### Testing correlation of 0.1 with 
## "arthritis", idiopathic", "inflammation", "inflammatory", 
## "juvenile", "polyarticular", "rheumatic", "rheumatoid" and "systemic"
findAssocs(dtmr, "idiopathic", 0.1)
findAssocs(dtmr, "inflammatory", 0.1)
findAssocs(dtmr, "polyarticular", 0.1)
findAssocs(dtmr, "rheumatic", 0.1)
findAssocs(dtmr, "rheumatoid", 0.1)
findAssocs(dtmr, "systemic", 0.1)
findAssocs(dtmr, "patient", 0.1)

# Basic graphics

## Frequency histogram (frequency > 100)
wf <- data.frame(term=names(freqr),occurrences=freqr)
p <- ggplot(subset(wf, freqr>100), aes(term, occurrences))
p <- p + geom_bar(stat="identity")
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1)) +
               labs(x="Termos", y="Ocorrências")
p

## Correlation plot
#cor_2 <- cor(as.matrix(dtmr))
#corrplot(cor_2, method = "color")

## Wordcloud
### Setting the same seed each time ensures consistent look across clouds
set.seed(42)
### Limit words by specifying min frequency
wordcloud(names(freqr),freqr, min.freq=70, colors=brewer.pal(6,"Dark2"))

# Cluster analysis
## Categorize documents into groups based on likeness (euclidean distance)  

## Hierarchical clustering
### Agglomerative: each document in its cluster and clusters are merged with the closest ones
### Divisive: each document in a single cluster and the cluster is split 

### Transform dtmr into a matrix
m <- as.matrix(dtmr)
### Compute distance between document vectors
d <- dist(m, method = "cosine")
### Run hierarchical clustering using Ward's method
groups <- hclust(d, method = "ward.D")
### Plot dendrogram; hang ensures labels fall below tree
plot(groups, hang = -1)
#### Too big :(
### Show clusters
#rect.hclust(groups, 6)

### Use WGCNA to try to automatically determine clusters based on dendrogram cutting
#dynamicMods <- cutreeDynamic(dendro=groups, minClusterSize = 200, cutHeight = 40)
dynamicMods <- cutree(tree=groups, h=28)
table(dynamicMods)
### Plot dendrogram and modules
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8, 6)
plotDendroAndColors(groups, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = -1,
                    addGuide = TRUE, main = "Dendrogram and module colors")

### Get module documents
modColors <- cbind(dtmr$dimnames$Docs, dynamicMods, dynamicColors)
modColors <- modColors[, -1]
#write.table(modColors, "docModules.txt", sep = "\t", quote=FALSE, row.names=F)

### Get wordclouds for each module
dtmrMat <- as.matrix(dtmr)
colors <- unique(modColors[, "dynamicColors"])
proper <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

#### Create labels file
globalLabel <- data.frame(docs=character(), color=character())

#pdf(file = "wordclouds.pdf",onefile=TRUE)
for(color in colors){
    print(color)
    colorMat <- subset(modColors, subset = dynamicColors == color)
    dtmrColor <- dtmr[rownames(colorMat), ]
    
    # Create label file
    colorLabel <- cbind(rownames(dtmrColor), color)
    globalLabel <- rbind(globalLabel, colorLabel)
    
    freqrMat <- colSums(as.matrix(dtmrColor))
    length(freqrMat)
    ## Create sort order (desc)
    ordrMat <- order(freqrMat, decreasing=TRUE)
    ## Inspect most and least frequently occurring terms
    print(freqrMat[head(ordrMat)])
    print(freqrMat[tail(ordrMat)])
    
    #wf <- data.frame(term=names(freqrMat),occurrences=freqrMat)
    #p <- ggplot(subset(wf, freqrMat>100), aes(term, occurrences))
    #p <- p + geom_bar(stat="identity")
    #p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
    #p
    
    ## Wordcloud
    ### Setting the same seed each time ensures consistent look across clouds
    set.seed(42)
    ### Limit words by specifying min frequency
    #layout(matrix(c(1, 2), nrow=3), heights=c(1, 4))
    par(mar=rep(0, 4))
    plot.new()
    #text(x=0.5, y=0.5, proper(color))
    pdf(paste0("wordcloud_", color, ".pdf"))
    wordcloud(names(freqrMat), main=color,freqrMat, min.freq=70, colors=brewer.pal(6,"Dark2"))
    dev.off()
    
    ## Word correlation plots
    #par(par.old)
    corColor <- cor(as.matrix(dtmrColor))
    p = corColor
    p[is.na(corColor)] <- 0.2
    p[is.na(corColor)==F] <- 0
    corColor[is.na(corColor)] <- 0
    #corrplot(corColor, method="color", is.corr=T, p.mat=p, sig.level=0.1, order="FPC")
    
    # Create adjacency matrices
    mColor <- as.matrix(dtmrColor)
    dColor <- dist(mColor, method = "cosine")
    n <- 75
    percent_n <- n*(length(as.vector(dColor))/100)
    element_n <- sort(dColor)[percent_n]
    dColor[dColor < element_n] <- 0
    dColor <- round(dColor, 3)
    #write.table(as.matrix(dColor),file=paste0(color, "_adjacencyMatrix.csv"), sep=";", quote=FALSE, col.names=NA)
    
}
#dev.off()
names(globalLabel) <- c("docs", "color")
write.table(globalLabel, "globalLabel.csv", sep=";", row.names=FALSE)

## PCA
m2 <- cbind(m, modColors[, "dynamicColors"])
colnames(m2)[length(colnames(m2))] <- "dynamicColors"
x <- as.data.frame(m2)

autoplot(prcomp(m), data=m2, colour="dynamicColors") + 
    scale_color_manual(values=levels(x$dynamicColors)) +
    theme(panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank()) +
    labs(title="PCA")

## K-means (takes too long and isn't good)
### 2 clusters, 100 starting configurations
#kfit <- kmeans(d, 2, nstart=100)
#clusplot(m, kfit$cluster, color=T, shade=T, labels=2, lines=0, col.p=kfit$cluster)

### Try to determine the best number of clusters (elbow method)
#### Look for 'elbow' in plot of summed intra-cluster distances (withinSS) as fn of k 
#### TAKES TOO LONG!!
#wss <- 2:20
#for (i in 2:20) wss[i] <- sum(kmeans(d, centers=i, nstart=25)$withinss)

# LDA
## Try to infer the latent topic structure given the words and document
## Recreate the documents in the corpus by adjusting the relative importance 
## of topics in documents and words in topics iteratively.

##Set parameters for Gibbs sampling
#l <- list()
#l[c("burnin", "iter", "thin", "seed", "nstart", "best", "k")] <- list(4000, 2000, 500, list(2003,5,63,100001,765), 5, TRUE, 5)
#list2env(setNames(l,names(l)), envir = parent.frame()) 
### Number of topics 

## Run LDA using Gibbs sampling (first, on the whole corpus!)
### Doesn't work
#ldaOut <-LDA(dtmr,k, method="Gibbs", 
#             control=list(nstart=nstart, seed = seed, best=best, 
#                          burnin = burnin, iter = iter, thin=thin))



# Visualizing document similarities as networks
## We already computed the cosine similarities between documents using the dist function
## Convert distance matrix diagonal elements to 0
n <- 90
dDist <- 1 - as.matrix(d)
percent_n <- n*(length(as.vector(dDist))/100)
element_n <- sort(as.vector(dDist))[percent_n]
dDist[dDist < element_n] <- 0
dDist <- round(dDist, 3)
diag(dDist) <- 0
#write.table(as.matrix(dDist),file="AdjacencyMatrix.txt", sep=";", quote=FALSE, col.names=NA)

## Get graph edgelist
g <- graph.adjacency(as.matrix(dDist), weighted=TRUE)
all_edges <- get.data.frame(g)
dim(all_edges)
names(all_edges) <- c("Source", "Target", "Weight")
write.table(all_edges,file="all_edges.csv", sep=";", quote=FALSE, row.names=FALSE)

# Redo wordcloud only with most important words
imp_words <- read.delim("best_class_words.txt", stringsAsFactors = FALSE)
imp_words <- imp_words[,1]

## Subset dtmr with the important words only (re-ran dtmr WITHOUT STEMMING!!!!!)
imp_dtmr <- dtmr[,intersect(colnames(dtmr), imp_words)]

## Calculate frequencies
imp_freq <- colSums(as.matrix(imp_dtmr))

## Plot wordcloud
set.seed(42)
### Limit words by specifying min frequency
png("wordcloud_imp.png")
wordcloud(names(imp_freq), imp_freq, min.freq=70, colors=brewer.pal(6,"Dark2"))
dev.off()

