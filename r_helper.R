

# https://rdrr.io/rforge/Vennerable/man/Venn.html#heading-2


library(Vennerable)

# read table
df <- read.table('enhancerlists.txt', header = TRUE)

# get sample names
samples <- colnames(df)

# find randoms
# be careful, this is hardcoded length vector. if the number of random samples
# changes, this must as well
subVennVec = rep(0,32)
print(subVennVec)
for (n in 1:10){
  print(n)
  subs <- sample(samples,5)
  subVenn <- Venn(df[,subs])
  subVennVec= subVennVec+c(Weights(subVenn))

}

# get average
subVennVec <- round(subVennVec/10)

# create new Venn object
meta <- Venn(SetNames=c('A','B','C','D','E'))
Weights(meta) <- subVennVec

# print meta
meta

# chow ruskey
pdf('output.pdf')
plot(meta,type='ChowRuskey')
dev.off()
