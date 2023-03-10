#' Corpus from raw document
#'
#' @param doc_raw the raw matrix of the documents
#' @param K number of topic
#'
#' @return list of corpus and vocabulary
#' @export
Corpus_alternative_word <- function(doc_raw, K){

  f <- lapply(doc_raw, function(d){
    dd <- tolower(d)
    dd <- tm::removePunctuation(dd)
    dd <- tm::removeNumbers(dd)
    dd <- unlist(strsplit(dd, " "))
    dd <- dd[!(dd %in% tm::stopwords())]
    dd <- tm::stemDocument(dd)
    dd <- tm::stripWhitespace(dd)
    dd <- stringr::str_subset(dd, ".+")
    return(dd)
  })

  vocabulary <- unique(unlist(lapply(f, function(x) unique(x))))
  V <- length(vocabulary)
  vocabulary <- data.frame(ID=1:V,word=vocabulary)

  corpus <- unlist(f)

  # length of the corpus:
  w <- numeric(length(corpus))
  # Total number of words:
  N <- length(corpus)
  # Number of documents:
  D <- length(doc_raw)

  # w is a vector length as the number of words in corpurs
  # for each word we have an index that connects it to the word in the vocabulary
  for(i in 1:N){
    w[i] <- vocabulary[which(corpus[i] == vocabulary[,2]),1]
  }

  # length of each document:
  length_docs <- unlist(lapply(f, function(ff){
    length(ff)
  }))
  sum((length_docs)) == N

  doc <- NULL
  for(d in 1:D){
    doc <- c(doc, rep(d, length_docs[d]))
  }

  Corpus_alternative <-  data.frame("Word"=w, "Doc"=doc)
  return(list("corpus"=Corpus_alternative, "vocabulary"=vocabulary))

}
