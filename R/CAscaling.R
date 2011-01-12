CAscaling <- function(Data, Range=3){
 nitem   <- ncol(Data)
 ncat    <- max(Data) + 1
 npers   <- nrow(Data)
 Combn   <- combn(nitem,2)
 Combn2  <- cbind(combn(nitem,2),rbind(Combn[2,],Combn[1,]))
 npairs  <- ncol(Combn)
 nitem   <- ncol(Data)
 Results <- matrix(0,4,nitem)
 
 Step1 <- matrix(0,nitem,ncol(Combn2))
 for(j in 1:nitem)
  for(k in unique(Data[,j])) 
   if(length(which(Data[,j] == k)) > 1){
    Cov <- cov(Data[which(Data[,j] == k),])
    ItemPairs <- which(Cov[row(Cov) > col(Cov)] < 0)
    if(length(ItemPairs) > 0)
     for(i in 1:length(ItemPairs)){
      Step1[Combn[2,ItemPairs[i]],which(Combn2[1,] == Combn[1,ItemPairs[i]] & Combn2[2,] == j)] <- 
        Step1[Combn[2,ItemPairs[i]],which(Combn2[1,] == Combn[1,ItemPairs[i]] & Combn2[2,] == j)] + length(which(Data[,j] == k))/nrow(Data)
      Step1[Combn[1,ItemPairs[i]],which(Combn2[1,] == Combn[2,ItemPairs[i]] & Combn2[2,] == j)] <- 
        Step1[Combn[1,ItemPairs[i]],which(Combn2[1,] == Combn[2,ItemPairs[i]] & Combn2[2,] == j)] + length(which(Data[,j] == k))/nrow(Data)
    }
  }
  
  Boxplot1 <- apply(Step1,2,sum)
  
  UpperFence1 <- max(Range * (quantile(Boxplot1,0.75) - quantile(Boxplot1,0.5)) + quantile(Boxplot1,0.75),(0.025*nitem) + quantile(Boxplot1,0.75))
  Remove1 <- unique(Combn2[1,which(Boxplot1 > UpperFence1)])
  Results[1,Remove1] <- Results[1,Remove1] + 1 
  

 S       <- var(Data)
 diag(S) <- 0
 Smax <- var(apply(Data, 2, sort))
 diag(Smax) <- 0
 Restcov <- order(apply(S, 1, sum)/apply(Smax, 1, sum))

 
 R.mat <- matrix(1, nitem, npairs)
 for(i in 1:npairs)
  R.mat[Combn[1:2,i],i] <- 0
 R <- as.matrix(Data) %*% R.mat
 
 Step2 <- matrix(0,nitem,nitem)
 Step3 <- rep(0,ncol(Combn))
 for(i in 1:npairs)
  for(j in unique(R[,i])) 
   if(length(which(R[,i] == j)) > 1)
    if(cov(Data[which(R[,i] == j),Combn[1,i]], Data[which(R[,i] == j),Combn[2,i]]) < 0){
     Step2[Combn[1,i],Combn[2,i]] <- Step2[Combn[1,i],Combn[2,i]] + length(which(R[,i]  == j))/nrow(Data)
      Step2[Combn[2,i],Combn[1,i]] <- Step2[Combn[2,i],Combn[1,i]] + length(which(R[,i]  == j))/nrow(Data)
     Step3[which(Combn[1,] == Combn[1,i] & Combn[2,] == Combn[2,i])] <- 
      Step3[which(Combn[1,] == Combn[1,i] & Combn[2,] == Combn[2,i])] + length(which(R[,i]  == j))/nrow(Data)
    }   

     
   Boxplot2 <- apply(Step2,2,sum)
   UpperFence2 <- max(Range * (quantile(Boxplot2,0.75) - quantile(Boxplot2,0.5)) + quantile(Boxplot2,0.75),(0.025*nitem) + quantile(Boxplot2,0.75))
   Remove2 <- which(Boxplot2 > UpperFence2)
   Results[2,Remove2] <- Results[2,Remove2] + 1 

  
  Boxplot3 <- Step3
  UpperFence3 <- max(Range * (quantile(Boxplot3,0.75) - quantile(Boxplot3,0.5)) + quantile(Boxplot3,0.75),(0.0025*nitem) + quantile(Boxplot3,0.75))
  Remove3 <- as.matrix(Combn[,which(Boxplot3 > UpperFence3)])

  if(length(Remove3) > 0)
   repeat{
    TableRem <- table(Remove3)
    if(length(Remove3) == 0) break
    if(max(TableRem) > 1){
     Items <- as.numeric(names(which(TableRem == max(TableRem))))
     Results[3,Items] <- Results[3,Items] + 1
     for(i in 1:length(Items))
      Remove3 <- as.matrix(Remove3[,which(apply(Remove3,2,function(x) any(x==Items[i])) == F)])
     next
    }
    if(max(TableRem) == 1){
     for(i in 1:ncol(Remove3)){
      HItem <- rep(0,2)
      HItem[1] <- which(Restcov == Remove3[1,i])
      HItem[2] <- which(Restcov == Remove3[2,i])
      Results[3,Remove3[which(HItem==min(HItem)),i]] <- Results[3,Remove3[which(HItem==min(HItem)),i]] + 1
     }
     break
    }
   }
  
  
 Results[4,] <- apply(Results[1:3,],2,sum)
 
 return(list(Output = Results,Scores = list(Boxplot1,Boxplot2,Boxplot3),Fences = list(UpperFence1,UpperFence2,UpperFence3)))

}
