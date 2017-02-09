
# info on vennerable usage
# https://rdrr.io/rforge/Vennerable/man/Venn.html#heading-2

library(Vennerable)

# get cmd line args
args <- commandArgs(TRUE)
sampleSize <- as.numeric(args[1])
vennVecLen <- 2^sampleSize
iters <- as.numeric(args[2])
samples <- unlist(strsplit(args[4],','))
cats <- unlist(strsplit(args[5],','))

randoms <- vector("list",iters)
for (n in 1:iters){
  print('RANDOMIZING')
  subs <- sample(samples,sampleSize)
  randoms[[n]] <- subs
}


plotter <- function(rnds,ctg) {
  print("plotting to:")
  print(paste(ctg,'_output.pdf',sep=''))
  # read table
  df <- read.table(paste(ctg,'_enhancerlists.txt',sep=''), header = TRUE)

  subVennVec = rep(0,vennVecLen)
  for (n in 1:iters){
    print("iterating")
    # get the names of samples as they appear in table
    real_names <- colnames(df)[grep(paste(rnds[[n]],collapse="|"),colnames(df))]
    subVenn <- Venn(df[,real_names])
    subVennVec= subVennVec+c(Weights(subVenn))
  }

  # get average
  subVennVec <- round(subVennVec/iters)
  covered_enhancers <- sum(subVennVec)
  subVennVec <- round((subVennVec/covered_enhancers*100,2)\

  # create new Venn object with just letters
  meta <- Venn(SetNames=letters[1:sampleSize])

  # populate it with the averaged percentage weights
  Weights(meta) <- subVennVec

  # print meta
  meta

  # chow ruskey
  pdf(paste(ctg,'_output.pdf',sep=''))
  plot(meta,type='ChowRuskey')
  dev.off()

}

for (c in cats) {
  plotter(randoms,c)
}
