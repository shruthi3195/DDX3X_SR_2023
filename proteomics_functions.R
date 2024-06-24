
#Quick data functions:
pD_channel <- function(df){
  
  df <- df %>% dplyr::mutate("channel_name" = ifelse(grepl("-0|mTRAQ0", Modified.Sequence), "mTRAQ0",
                                                     ifelse(grepl("-4|mTRAQ4", Modified.Sequence), "mTRAQ4", "mTRAQ8")))
  
  return(df)
}

pD_seqcharge <- function(df){
  df$seqcharge <- paste0(df$Modified.Sequence, df$Precursor.Charge)
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-8\\)")
  return(df)
}


# K - nearest neighbors imputation

hknn<-function(dat, k){
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))
  
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      #print(length(closest.columns))
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
        
      }
      
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ] )
        
      }
      
      
    }
    
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
    
  }
  
  return(dat.imp)
  
}