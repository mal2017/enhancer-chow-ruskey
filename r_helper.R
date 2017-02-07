
# info on vennerable usage
# https://rdrr.io/rforge/Vennerable/man/Venn.html#heading-2

library(Vennerable)

# get cmd line args
args <- commandArgs(TRUE)
sampleSize <- as.numeric(args[1])
vennVecLen <- 2^sampleSize
iters <- as.numeric(args[2])

# read table
df <- read.table('enhancerlists.txt', header = TRUE)

# get sample names
samples <- colnames(df)

# find randoms
#
subVennVec = rep(0,vennVecLen)
print(subVennVec)
for (n in 1:iters){
  print(n)
  subs <- sample(samples,sampleSize)
  subVenn <- Venn(df[,subs])
  subVennVec= subVennVec+c(Weights(subVenn))

}

# get average
subVennVec <- round(subVennVec/iters)

# create new Venn object
meta <- Venn(SetNames=letters[1:sampleSize])
Weights(meta) <- subVennVec

# print meta
meta

# chow ruskey
pdf(args[3]+'/output.pdf')
plot(meta,type='ChowRuskey')
dev.off()
