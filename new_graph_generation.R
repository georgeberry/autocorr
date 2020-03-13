###Create simulation of network that has N, nodes/20000 edges.

# Algorithm
# Generate nodes
#   - 
# 


###Population Attributes
#Size of Classes - Size of X
#Distribution of Y and Z
#Relationship between X and Z.
#Relationship between Y and Z.

N = 100
E = 500
LatentNoise = 0.5
###Things to vary:
#Degree distribution/Preferential Attachment
#Homophily along X (Structural Homophily)
#Homophily along Y (Latent Homophily)

########
PrefAttach_Mod = 0
Struct_Homoph_Mod = 1
Latent_Homoph_Mod = 0

###Node attributes:
#X Class
#Y Latent Variance
#Z Observed Variance

###For Loops That Generate X,Y,Z
NodeMatrix <- matrix(nrow=0,ncol=5)
for(a in 1:N){
  x <- sample(0:1,1)
  y <- rnorm(1,0,LatentNoise) + x 
  z <- rnorm(1,0,1) + x
  NodeEdges <- 0
  NodeRow <- cbind(a,x,y,z,NodeEdges)
  #print(NodeRow)
  NodeMatrix <- rbind(NodeMatrix,NodeRow)
}
NodeMatrix<-as.data.frame(NodeMatrix)
###Visualize Difference in Node Class and Predictive Value
#plot(density(NodeMatrix[which(NodeMatrix$x==0),]$y),col='red')
#lines(density(NodeMatrix[which(NodeMatrix$x==1),]$y),col='blue')
#Classify for X based on Y Using Simple Logistic Regression
Prediction <-as.numeric(predict(glm(formula = x ~ y, family = "binomial", data = NodeMatrix))>0.5)
NodeMatrix<-cbind(NodeMatrix,Prediction)
#Generate Edges
#Need to Multiple Edges, Loops
EdgeList <-matrix(nrow=0,ncol=8)
while(nrow(EdgeList) < E){
  ProposedEdge <- NodeMatrix[sample(nrow(NodeMatrix), 2),]
  #print(ProposedEdge)
  Node_1 =ProposedEdge[1,1]
  Node_2 =ProposedEdge[1,2]
  X_1 <-ProposedEdge[1,2]
  X_2 <- ProposedEdge[2,2]
  Y_1 <- ProposedEdge[1,3]
  Y_2 <- ProposedEdge[2,3]
  Pred_1 <- ProposedEdge[1,6]
  Pred_2 <- ProposedEdge[2,6]
  X_diff <- abs(X_1-X_2)
  Y_diff <- abs(Y_1-Y_2)
  X_edges <- ProposedEdge[1,5]
  Y_edges <- ProposedEdge[2,5]
  #print(c(X_diff,Y_diff,X_edges,Y_edges))
  logit = (X_edges + Y_edges)*PrefAttach_Mod - Struct_Homoph_Mod*X_diff - Latent_Homoph_Mod*Y_diff
  EdgeProb = logit/(1+logit)
  RandomNumber = runif(1,0,1)
  #print(RandomNumber,EdgeProb)
  if(RandomNumber > EdgeProb){
    NewEdge = cbind(Node_1,Node_2,X_1,Y_1,X_2,Y_2,Pred_1,Pred_2)
    EdgeList= rbind(EdgeList,NewEdge)
    NodeMatrix$NodeEdges[Node_1] = X_edges+1
    NodeMatrix$NodeEdges[Node_2] = Y_edges+1
  }
}
EdgeList<-as.data.frame(EdgeList)

ClassificiationTable = table(NodeMatrix$x,NodeMatrix$Prediction)
TrueHomophilyTable = table(EdgeList$X_1,EdgeList$X_2)
ClassifiedHomphilyTable = table(EdgeList$Pred_1,EdgeList$Pred_2)
print(ClassificiationTable)
print(TrueHomophilyTable)
print(ClassifiedHomphilyTable)

ClassifierAccuracy = (ClassificiationTable[1,1]+ClassificiationTable[2,2])/N
ActualHomophily = (TrueHomophilyTable[1,1]+TrueHomophilyTable[2,2])/E
ClassifiedHomophily = (ClassifiedHomphilyTable[1,1]+ClassifiedHomphilyTable[2,2])/E
print(c(ClassifierAccuracy,ActualHomophily,ClassifiedHomophily))



N = 100
E = 500
LatentNoise = 0.5
###Things to vary:
#Degree distribution/Preferential Attachment
#Homophily along X (Structural Homophily)
#Homophily along Y (Latent Homophily)

########
PrefAttach_Mod = 0
Struct_Homoph_Mod = 1
Latent_Homoph_Mod = 0

Simulate_Homophily <- function(
  N,
  E,
  LatentNoise,
  PrefAttach_Mod,
  Struct_Homoph_Mod,
  Latent_Homoph_Mod,
  Iterations
){
  FinalMatrix = matrix(nrow=0,ncol=4)
  for(c in 1:Iterations){
    NodeMatrix <- matrix(nrow=0,ncol=5)

    # this generates nodes charaacteristics, x y z
    for(a in 1:N){
      x <- sample(0:1,1) # samples a 1/0 X
      y <- rnorm(1,0,LatentNoise) + x 
      z <- rnorm(1,0,1) + x
      NodeEdges <- 0
      NodeRow <- cbind(a,x,y,z,NodeEdges)
      #print(NodeRow)
      NodeMatrix <- rbind(NodeMatrix,NodeRow)
    }
    NodeMatrix<-as.data.frame(NodeMatrix)
    ###Visualize Difference in Node Class and Predictive Value
    #plot(density(NodeMatrix[which(NodeMatrix$x==0),]$y),col='red')
    #lines(density(NodeMatrix[which(NodeMatrix$x==1),]$y),col='blue')


    # predictions

    #Classify for X based on Y Using Simple Logistic Regression
    Prediction <-as.numeric(predict(glm(formula = x ~ y, family = "binomial", data = NodeMatrix))>0.5)
    NodeMatrix<-cbind(NodeMatrix,Prediction)
    #Generate Edges


    # now add edges

    #Need to Multiple Edges, Loops
    EdgeList <-matrix(nrow=0,ncol=8)
    while(nrow(EdgeList) < E){

      # sample a node
      ProposedEdge <- NodeMatrix[sample(nrow(NodeMatrix), 2),]


      #print(ProposedEdge)
      Node_1 =ProposedEdge[1,1]
      Node_2 =ProposedEdge[1,2]

      X_1 <-ProposedEdge[1,2]
      X_2 <- ProposedEdge[2,2]
      Y_1 <- ProposedEdge[1,3]
      Y_2 <- ProposedEdge[2,3]
      Pred_1 <- ProposedEdge[1,6]
      Pred_2 <- ProposedEdge[2,6]

      # what is this? pairwise difference in x and y?
      X_diff <- abs(X_1-X_2)
      Y_diff <- abs(Y_1-Y_2)

      X_edges <- ProposedEdge[1,5]
      Y_edges <- ProposedEdge[2,5]
      #print(c(X_diff,Y_diff,X_edges,Y_edges))
      logit = (
        1 +
        (X_edges + Y_edges) * PrefAttach_Mod -
        Struct_Homoph_Mod * X_diff -
        Latent_Homoph_Mod * Y_diff
      )
      EdgeProb = logit/(1+logit)
      RandomNumber = runif(1,0,1)
      #print(RandomNumber,EdgeProb)
      if(RandomNumber > EdgeProb){
        NewEdge = cbind(Node_1,Node_2,X_1,Y_1,X_2,Y_2,Pred_1,Pred_2)
        EdgeList= rbind(EdgeList,NewEdge)
        NodeMatrix$NodeEdges[Node_1] = X_edges+1
        NodeMatrix$NodeEdges[Node_2] = Y_edges+1
      }
    }
    EdgeList<-as.data.frame(EdgeList)
    ClassificiationTable = table(NodeMatrix$x,NodeMatrix$Prediction)
    TrueHomophilyTable = table(EdgeList$X_1,EdgeList$X_2)
    ClassifiedHomphilyTable = table(EdgeList$Pred_1,EdgeList$Pred_2)
    ClassifierAccuracy = (ClassificiationTable[1,1]+ClassificiationTable[2,2])/N
    ActualHomophily = (TrueHomophilyTable[1,1]+TrueHomophilyTable[2,2])/E
    ClassifiedHomophily = (ClassifiedHomphilyTable[1,1]+ClassifiedHomphilyTable[2,2])/E
    NewFinalRow <- c(c,ClassifierAccuracy,ActualHomophily,ClassifiedHomophily)
    FinalMatrix <-rbind(FinalMatrix,NewFinalRow)
  }
  FinalMatrix<-as.data.frame(FinalMatrix)
  names(FinalMatrix)<-c('Trial','Accuracy','True_H','Classified_H')
  return(FinalMatrix)
}

SampleResult1 <- Simulate_Homophily(200,1000,0.5,0,1,0,50)
SampleResult2 <- Simulate_Homophily(200,1000,0.5,0,0,1,50)

plot(SampleResult1$True_H,SampleResult1$Classified_H,col='Blue',xlim=c(0.4,0.7),ylim=c(0.4,0.7))
points(SampleResult2$True_H,SampleResult2$Classified_H,col='Red')

