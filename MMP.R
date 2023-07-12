# Modified maximum parsimony
#

#' Modified maximum parsimony - MMP
#'
#' @param file The file name. Need to be an .xlsx file. It is the output file of DEVOLUTION.
#' @param tumorname The name for your sample. This name will be used in the filename for the output file.
#'
#' @return The phylo object for the MMP-tree.
#' @export
#'
#' @examples MMP(filename="DEVOLUTION.xlsx",tumorname="T1")
MMP <- function(file,tumorname){

  load_matrix <- function(filename, sheetname) {
    data <- as.data.frame(read_xlsx(filename, sheetname)) #Reading the xlsx file and saving it in the variable data.
    subdata <- data[ c(1:nrow(data)), c(1:ncol(data)) ] #Extracting the part of the file we are interested in.
    subdata <- subdata[is.na(subdata[,1])==FALSE,]
    return(subdata)
  }

  EM_full <- load_matrix(filename=file,sheetname="Event matrix")
  if(colnames(EM_full)[ncol(EM_full)]!="Type"){ #If the type column is missing we add a default that all are whole chromosome alterations.
    typecol <- as.matrix(rep("W",nrow(EM_full)))
    colnames(typecol) <- "Type"
    EM_full <- cbind(EM_full,typecol)
  }
  EM <- EM_full[,1:(ncol(EM_full)-1)]
  type <- EM_full[,ncol(EM_full)] #This is transformed compared with the one from DEV.
  pies <- load_matrix(filename=file,sheetname="Pies")
  pies <- rbind(colnames(pies),pies)
  overview <- load_matrix(filename=file,sheetname="Overview")
  overview <- overview[,2:ncol(overview)]

  #Create the tip labels.
  tip.label <- c(colnames(EM)[2:length(colnames(EM))])

  #The edge matrix.
  edge <- matrix(0,length(tip.label)*4,2)

  #The edge length matrix.
  edge.length <- t(matrix(0,length(tip.label)*4,1))

  #Base matrix. The subclones allocated at the same height as the node.
  base <- t(matrix(0,length(tip.label),2))

  #Ordering the EM after the rowsum.
  EM_or <- order(rowSums(EM[,2:ncol(EM)]),decreasing=TRUE)
  EM <- EM[EM_or,]

  #Identifying alterations involved in or responsible for contradictions in the tree.
  i <- 2
  M <- list()
  contr <- vector()
  cause <- vector()
  m <- 1
  for(i in 2:nrow(EM)){ #Looping over the events (rows) in the EM.
    diff1 <- EM[i-1,2:ncol(EM)]
    diff2 <- EM[i,2:ncol(EM)]
    diff3 <- length(diff1[diff1 != diff2]) #Checking if the two events differ from one another.

    if(as.numeric(diff3)!=0){ #If they do not differ from one another, they do not cause contradictions.
      ev_m <- matrix(0,nrow(EM),ncol(EM)-1) #Here we will save the mothers of each alteration.

      j <- 1
      for(j in 2:ncol(EM)){ #Looping over the presence of this alteration across subclones.

        if(EM[i,j]=="1"){
          #The alteration is present in this subclone.
          EM_sub <- EM[1:(i-1),] #Extracting the rows above the alteration in the matrix, i.e. the alterations in which it is nested.
          EM_sub <- EM_sub[EM_sub[,j]!="0",c(1,j)] #Extracting a particular column.
          ev_m[1:nrow(EM_sub),j] <- EM_sub[,1] #Saving the alterations present above the alteration in this column.
        }

        j <- j+1
      }

      #Same mother in all cases?
      ev_m <- ev_m[,ev_m[1,]!="0"] #Removing the empty columns.

      if(is.null(dim(ev_m))==FALSE){
        #eq <- all(ev_m%in%ev_m[,1]) #Checking if all columns are equal.
        eq <- all(ev_m==ev_m[,1]) #Checking if all columns are equal. #Maybe this is wrong. They do not have to have the same position right?
        # print(ev_m)
        # print(eq)
        if(eq==FALSE){
          #The allocation differs across samples.
          #Checking if the events that differ between samples are of equal size to the alteration in question.

          t <- table(ev_m) #Looking at in how many of the columns the events are present.
          u <- t[t<ncol(ev_m)] #Identifying alterations that are "new"/unique across the columns. These are involved in the contradiction.

          #Checking if it is equal in size to the daughtersubclone.
          pos_m <- which(overview[,1]%in%names(u))
          pos_d <- which(overview[,1]%in%EM[i,1])

          k <- 1
          for(k in 1:length(pos_m)){
            diff_md <- as.numeric(overview[pos_m[k],2:ncol(overview)])-as.numeric(overview[pos_d,2:ncol(overview)]) #Comparing their sizes across samples.
            # print(EM[i,1])
            # print(diff_md)

            neg <- which(diff_md<0)
            posi <- which(diff_md>0)

            if(length(neg)>0&&length(posi)>0){
              print("The events are crossing each other")
              print(names(u)[k])
              print(EM[i,1]) #this one is causing it.

              print("These events are also involved")
              inv <- EM[(i+1):nrow(EM),c(1,which(EM[3,]==1))]
              print(inv[,1])

              cause <- c(cause,EM[i,1])
              contr <- c(contr,names(u[k]),EM[i,1],inv[,1])

            }

            k <- k+1
          }

        }
      }

      lst <- list(events_mother = ev_m, event_name = EM[i,1])
      M[[m]] <- lst #Saving the matrix with the mothers.
      m <- m+1

    }
    i <- i+1
  }
  contr <- unique(contr)
  cause <- unique(cause)
  print("This event will cause contradiction in the tree. Revise?")
  print(cause)
  print("These events will also be involved due to this.")
  print(contr)


  i <- 2
  s <- 1
  d <- 1
  lv <- 1
  diff_m <- list()
  for(i in 2:ncol(EM)){ #Looping over the subclones.
    print(i)
    name <- tip.label[1]
    if(i==2){ #Allocating the stem and the normal cell.
      #The stem
      nodenr <- length(tip.label)+1
      edge[s,1] <- nodenr
      edge[s,2] <- i-1
      edge.length[s] <- 0 #The stem.
      allocated <- EM[,i]
      remaining <- EM[,-i] #Removing it.
      s <- s+1

      #The Normal cell.
      edge[s,1] <- length(tip.label)+1
      edge[s,2] <- which(tip.label=="Normal")
      edge.length[s] <- length(EM[EM[,2]=="1",2]) #The length of the stem.
      allocated <- cbind(allocated,remaining[,ncol(remaining)])
      colnames(allocated)[1] <- colnames(EM)[i]
      colnames(allocated)[2] <- colnames(remaining)[ncol(remaining)]
      rownames(allocated) <- EM[,1]
      remaining <- remaining[,-ncol(remaining)] #Removing it.
      s <-s+1

    }else if(i>2 && i<ncol(EM)){ #Allocating the other subclones.

      j <- 1
      print(colnames(remaining))
      print(colnames(allocated))

      for(j in 1:ncol(allocated)){

        #Computing the diff matrix.
        if(ncol(remaining)!=2){
          diff <- sweep(remaining[,2:ncol(remaining)],1,allocated[,j])
        }else{
          diff <- as.matrix(remaining[,2:ncol(remaining)]-allocated[,j]) #If it is 2 we cannot use sweep.
          colnames(diff)[1] <- colnames(remaining)[2]
        }

        #Saving the diffmatrix and information about which subclone we compared to in a list.
        lst <- list(diffmatrix=diff, allocated = allocated, name = colnames(allocated)[j])
        diff_m[[d]] <- lst
        d <- d+1

        if(j!=1){
          min_sum <- colSums(abs(diff))
          min_sum_sub2 <- t(as.matrix(min_sum[which(min_sum==min(min_sum))]))
          name <- colnames(min_sum_sub2)
          min_sum_sub2 <- rbind(min_sum_sub2,colnames(allocated)[j]) #Adding the mother.
          min_sum_sub2 <- rbind(min_sum_sub2,name) #Adding the daughter name.


          name_pos <- which(colnames(diff)%in%name)
          n <- 1
          for(n in 1:length(name)){
            if(n != 1){
              neg_1 <- length(diff[diff[,name_pos[n]]==-1,name_pos[n]]) #Back mutations.
              neg <- c(neg,neg_1)
            }else{
              neg <- length(diff[diff[,name_pos[n]]==-1,name_pos[n]])
            }
            n <- n+1
          }

          min_sum_sub2 <- rbind(min_sum_sub2,neg)
          min_sum_sub1 <- cbind(min_sum_sub1,min_sum_sub2) #Combining it with the previous one.
        }else{
          min_sum <- colSums(abs(diff))
          min_sum_sub1 <- t(as.matrix(min_sum[which(min_sum==min(min_sum))]))
          name <- colnames(min_sum_sub1)
          min_sum_sub1 <- rbind(min_sum_sub1,colnames(allocated)[j]) #Adding the mother.
          min_sum_sub1 <- rbind(min_sum_sub1,name) #Adding the daughter name.
          neg <- length(diff[diff[,1]==-1,1]) #Back mutations.
          min_sum_sub1 <- rbind(min_sum_sub1,neg) #Back mutation indication.

        }

        j <- j+1
      }

      #Finding the smallest one which does not imply back mutations.
      print(min_sum_sub1)

      sub <- min_sum_sub1[,min_sum_sub1[4,]=="0"]
      sub <- as.matrix(sub[,which(as.numeric(sub[1,])==min(as.numeric(sub[1,])))]) #The one we should allocate and where.


      if(is.na(ncol(sub))==TRUE){
        #We cannot choose anyone that does imply back mutations.
        sub <- as.matrix(min_sum_sub1[,which(as.numeric(min_sum_sub1[1,])==min(as.numeric(min_sum_sub1[1,])))]) #The one we should allocate and where.
      }
      print("The sub")
      print(sub) #Diff, mother, daughter.

      #Checking if one of the new alterations in the new daughter subclone is one giving contradictions.
      mcol <- which(colnames(EM)==sub[2,1]) #mother
      dcol <- which(colnames(EM)==sub[3,1]) #daughter
      diff_dm <- as.numeric(EM[,dcol])-as.numeric(EM[,mcol])
      ev_dm <- EM[which(diff_dm!=0),1]
      cause_ev <- ev_dm[ev_dm%in%cause]

      #Checking if it can be placed in the current mother based on pie sizes.
      #Also checking if there are other clones with larger pie size that could be allocated instead.
      #Also checking if there are other subclones more suitable to allocate.

      #Daughter subclone pie sizes.
      daughtername <- sub[3,1]
      mothername <- sub[2,1]
      pie_pos_d <- which(pies[,1] == daughtername) #The position in the pie matrix for the daughterclone.
      pie_d <- pies[pie_pos_d:(pie_pos_d+1),] #Removing the columns with zeros.
      pie_d <- pie_d[,pie_d[1,]!="0"]
      pie_d <- pie_d[,pie_d[2,]!="0"]

      #Checking if there is some other not allocated subclone that should be allocated instead.
      m <- 2
      for(m in 2:length(colnames(remaining))){

        if(colnames(remaining)[m]!=sub[3,1] && colnames(remaining)[m]!="Normal" && colnames(remaining)[m]!="Stem"){ #We do not need to compare it to itself.

          pie_pos_m <- which(pies[,1]==colnames(remaining)[m])
          pie_m <- pies[pie_pos_m:(pie_pos_m+1),] #Removing the columns with zeros.
          pie_m <- pie_m[,pie_m[1,]!="0"]
          pie_m <- pie_m[,pie_m[2,]!="0"]

          p <- 2
          for(p in 2:ncol(pie_m)){ #Looping over the pies in one of the remaining ones.
            s1 <- word(pie_m[1,p],1) #Sample
            s1_size <- pie_m[2,p] #The size in that sample for the other subclone.

            s2_size <- pie_d[2,match(s1,word(pie_d[1,],1))] #The size of the daughter we wanted to allocate in that sample.

            if(length(s2_size)!=0){ #Maybe the subclone is not even present in that sample.
              if(as.numeric(s1_size)+as.numeric(s2_size)>100 && as.numeric(s1_size)>as.numeric(s2_size)){
                #The daughter's size and the mother's size in that sample is larger than 100 % and should be nested in one another.
                #The other not allocated subclone is also larger in size.
                daughtername <- colnames(remaining)[m] #daughter
                print("Another daughter first")
                print(s1_size)
                print(s2_size)
                print(daughtername)
                print(mothername)

              }
            }

            p <- p+1
          }
        }
        m <- m+1
      }


      #Daughter subclone pie sizes.
      pie_pos_d <- which(pies[,1] == daughtername) #The position in the pie matrix for the daughterclone.
      pie_d <- pies[pie_pos_d:(pie_pos_d+1),] #Removing the columns with zeros.
      pie_d <- pie_d[,pie_d[1,]!="0"]
      pie_d <- pie_d[,pie_d[2,]!="0"]

      #Checking if the daughter need to be placed somewhere else due to pie sizes.
      m <- 1
      for(m in 1:length(colnames(allocated))){ #Looping over the already allocated clones.

        if(colnames(allocated)[m]!="Normal" && colnames(allocated)[m]!="Stem" && colnames(allocated)[m]!=mothername && substr(colnames(allocated)[m],start=1,stop=2)!="lv"){
          pie_pos_m <- which(pies[,1]==colnames(allocated)[m])
          pie_m <- pies[pie_pos_m:(pie_pos_m+1),] #Removing the columns with zeros.
          pie_m <- pie_m[,pie_m[1,]!="0"]
          pie_m <- pie_m[,pie_m[2,]!="0"]

          p <- 2
          for(p in 2:ncol(pie_m)){ #Looping over the pies in the mother clone.
            s1 <- word(pie_m[1,p],1) #Sample
            s1_size <- pie_m[2,p] #The size of the other mother for that sample.

            s2_size <- pie_d[2,match(s1,word(pie_d[1,],1))] #The size of the daughter for that sample.

            if(length(s2_size)!=0){
              if(as.numeric(s1_size)+as.numeric(s2_size)>100 && as.numeric(s1_size)>as.numeric(s2_size)){
                #The daughter's size and the mother's size in that sample is larger than 100 % and should be nested in one another.
                #The daughter is also smaller in size.

                #It may although be the case that the previously chosen mother is also nested in this one, in which
                #case we do not need to nest this one down here further down in the phylogeny.
                #Let's check if the new mother is an ancestor the the current one.
                #Remark! I'm only checking one step down here. We should check multiple steps actually...

                m_new <- which(mothername==tip.label)
                m_old <- which(colnames(allocated)[m]==tip.label)

                m_new_node <- edge[which(m_new==edge[,2]),1] #Node
                m_old_node <- edge[which(m_old==edge[,2]),1] #Node

                c1 <- which(paste(m_new_node,m_old_node)==paste(edge[,1],edge[,2]))
                c2 <- which(paste(m_old_node,m_new_node)==paste(edge[,1],edge[,2]))

                if(length(c1) == 0 && length(c2) == 0){
                  #They are not connected.
                  print("They are not connected") #Hence we will change the mother.

                  mothername <- colnames(allocated)[m] #New potential mother.
                  print("Daughter at other location")
                  print(s1_size)
                  print(s2_size)
                  print(daughtername)
                  print(mothername)
                }else{
                  print("They are connected") #We will not change the mother.
                }

              }
            }
            p < p+1
          }

        }
        m <- m+1
      }

      sub <- matrix(0,3,1)
      sub[1,1] <- sum(abs(as.numeric(EM[,which(colnames(EM)==mothername)])-as.numeric(EM[,which(colnames(EM)==daughtername)])))
      sub[2,1] <- mothername
      sub[3,1] <- daughtername

      print("The final sub chosen for i")
      print(sub)
      print(i)

      #Equal events.
      mcol <- which(colnames(allocated)==sub[2,1]) #mother
      dcol <- which(colnames(remaining)==sub[3,1]) #daughter

      m <- 1
      eq_m <- matrix(0,3,ncol(allocated))
      for(m in 1:ncol(allocated)){ #Making a matrix showing how many events the daughter has in common with the other ones.
        eq_m[1,m] <- colnames(allocated)[m] #Mothername.
        eq_m[2,m] <- sum(as.numeric((as.numeric(allocated[,m])+as.numeric(remaining[,dcol]))==2)) #Nr of events in common.
        eq_m[3,m] <- length(which(as.numeric(remaining[,dcol]-allocated[,m])==-1)) #Counting backmutations.
        m <- m+1
      }

      #Selecting the ones already allocated and finding which one the subclone has more in common to.
      #eq_new <- as.matrix(eq_m[,which(eq_m[3,]=="0")]) #Selecting the ones not implying back mutations.
      eq_new <- as.matrix(eq_m[,which(eq_m[1,]!="Normal")]) #Selecting the ones not implying placing it before the stem. It is not reasonable.

      print("Before")
      print(eq_m)
      print(eq_new)

      if(length(eq_new)!=0){
        #There are clones left in the matrix.
        if(sub[2,1]%in%eq_new[1,]==FALSE){
          #The mother implied backmutations and is no longer in the matrix.
          #We will keep the other clones and see if they are a better fit.
          print("The mother implies back mutations.")

        }else{
          #The mother did not imply back mutations. But let's also investigate if there is a better place for it to be placed.
          print("The mother does not imply back mutations.")
          #Selecting events which it has more in common to.
          eq_new <- as.matrix(eq_m[,which(as.numeric(eq_m[2,])>as.numeric(eq_m[2,which(eq_m[1,]==sub[2,1])]))])
        }
      }else{
        #There is no allocation that does not imply back mutations.
      }

      if(length(eq_new)!=0){
        print("SMALL")
        #Ordering the subclones in sub_part.
        eq_new <- as.matrix(eq_new[,order(as.numeric(eq_new[3,]),decreasing=FALSE)])
        print(eq_new)
        #The subclone is similar to a part of another subclone.
        sub_part <- as.matrix(eq_new[,which(as.numeric(eq_new[2,])==max(as.numeric(eq_new[2,])))]) #The clone chosen that we will create a new level to.
        print(sub_part)

        #We need to change the order of the "mother" and the "sub_part"-mother if the sub_part is lower than the mother in the tree.
        sw1 <- sub
        sw2 <- sub_part

        #Numbers
        mothernr <- which(tip.label==sub[2,1]) #Mother.
        daughternr <- which(tip.label==sub[3,1]) #Daughter

        #The mother.
        if(substr(sub[2,1],start=1,stop=2)=="lv"){
          #The mother is a level.
          mothernode <- lv_m[2,which(lv_m[1,]==sub[2,1])]
        }else{
          mothernr <- which(tip.label==sub[2,1]) #Mother.
          motherrow <- which(edge[,2] == mothernr)
          mothernode <- edge[motherrow,1] #The node in which the mother is allocated.
          lengths <- edge.length[motherrow]
          motherrow <- motherrow[which(lengths==0)] #Treating the case where the mother is at multiple positions.
          mothernode <- mothernode[which(lengths==0)] #The node in which the mother is allocated.
        }

        #The subclone it is similar to.
        if(substr(sub_part[1,1],start=1,stop=2)=="lv"){
          #The subclone it is similar to is a level.
          subnode <- lv_m[2,which(lv_m[1,]==sub_part[1,1])]
        }else{
          subnr <- which(tip.label==sub_part[1,1]) #Subclone that the daughter is part of.
          subrow <- which(edge[,2]==subnr) #The row at which the subclone the daughter is part of is.
          subnode <- edge[subrow,1] #The node in which the subclone the daughter is part of is allocated.
        }

        subnode_pos <- which(edge[,2]==subnode) #Where the node connects to its mother node.
        sublength <- edge.length[subnode_pos] #The length between the node of the already allocated subclone and its mother.

        nodenr <- nodenr+1 #The new node.

        #Changing
        #Short new segment.
        ev1 <- eq_m[2,which(eq_m[1,]==sub[2,1])] #The number of events it has in common to its mother.
        ev2 <- sub_part[2,1]#eq_m[2,which(eq_m[1,]==sub_part[2,1])] #sub_part[2,1] #The number of events it has in common to the other subclone.

        edge.length[s] <- abs(as.numeric(ev1)-as.numeric(ev2))
        print(edge.length[s])

        edge[s,1] <- edge[subnode_pos,1] #The node the subclone connects to below.      Previously: mothernode
        edge[s,2] <- nodenr
        s <- s+1

        #Change the mother node and branch lengths of the already allocated subclone.
        edge[subnode_pos,1] <- nodenr #Changing the mothernode to the new level.
        Lofmother <- sum(as.numeric(allocated[,which(colnames(allocated)==sub_part[1,1])]))
        #The old length of this branch minus how many unique events the mother has that the new subclone does not.

        edge.length[subnode_pos] <- Lofmother-as.numeric(ev2) #New length of this small segment to the subclone.

        #Allocating the new daughter clone.
        #Daughternode - Newnode
        edge[s,1] <- nodenr
        nodenr <- nodenr+1
        edge[s,2] <- nodenr

        Lofdaughter <- sum(as.numeric(remaining[,which(colnames(remaining)==sub[3,1])]))
        edge.length[s] <- Lofdaughter-as.numeric(sub_part[2,1]) #Length of the daughter minus what it had in common to the other mother.

        #Adding the new node level to the allocated matrix.
        print("Adding a new level")
        m1 <- as.numeric(allocated[,which(colnames(allocated)==sub_part[1,1])])
        d1 <- as.numeric(remaining[,which(colnames(remaining)==sub[3,1])])
        n1 <- as.matrix(cbind(m1,d1))

        allocated_vector_level <- as.matrix(as.numeric(rowSums(n1)==2))

        allocated <- cbind(allocated,allocated_vector_level) #The events a subclone at this level would have.
        nm <- paste("lv",nodenr-1,sep="") #The name.
        colnames(allocated)[ncol(allocated)] <- nm
        print(nm)

        #If we create a new level we also want to save the level name and the node belonging to each level.
        if(lv==1){
          #The first level we have created.
          lv_m <- as.matrix(c(nm,nodenr-1))
          lv <- lv+1
        }else{
          lv_m <- cbind(lv_m,as.matrix(c(nm,nodenr-1)))
        }

        s <- s+1
      }else{
        #It should be after this specific subclone.
        #Allocating it.
        mothernr <- which(tip.label==sub[2,1]) #Mother.
        daughternr <- which(tip.label==sub[3,1]) #Daughter

        motherrow <- which(edge[,2] == mothernr)
        mothernode <- edge[motherrow,1] #The node in which the mother is allocated.

        lengths <- edge.length[motherrow]
        motherrow <- motherrow[which(lengths==0)] #Treating the case where the mother is at multiple positions.
        mothernode <- mothernode[which(lengths==0)]

        #Daughternode - Mothernode
        edge[s,1] <- mothernode
        nodenr <- nodenr+1
        edge[s,2] <- nodenr

        edge.length[s] <- as.numeric(sub[1,1])
        s <- s+1

      }

      #Telling the user if we have instances of parallel evolution.
      daughtercol <- which(colnames(remaining)==sub[3,1])
      mothercol <- which(colnames(allocated)==tip.label[mothernr])
      diff_dm <- as.numeric(remaining[,daughtercol])-as.numeric(allocated[,mothercol]) #The unique alterations for this one.
      pos_u <- which(diff_dm==1) #The positions of the unique one. If it is already allocated at another place, it will be parallel.
      parallel <- which(allocated[pos_u,]==1,arr.ind = T) #Are the events present in other subclones allocated?

      #Daughter - Daughternode
      edge[s,1] <- nodenr
      edge[s,2] <- daughternr
      edge.length[s] <- 0
      s <- s+1

      #Change the allocated and remaining matrices.
      daughtercol <- which(colnames(remaining)==sub[3,1])
      allocated <- cbind(allocated,remaining[,daughtercol])
      colnames(allocated)[ncol(allocated)] <- sub[3,1]
      remaining <- remaining[,-daughtercol]

    }

    i <- i+1
  }


  #Tip.label
  tip.label

  #Edge.length
  edge.length <- edge.length[1:nrow(edge)]

  #Edge
  edge <- edge[edge[,1]!=0,] #Removing empty rows.
  edge.length <- edge.length[1:nrow(edge)] #Shortening the edge.length as well.


  #Edges should not have a node.
  nodes <- table(edge[,1])
  few <- as.numeric(names(nodes[which(nodes==1)])) #These nodes do only have one subclone.
  row_few <- which(edge[,1]%in%few)
  i <- 1
  for(i in 1:length(row_few)){

    row <- which(edge[,1]==few[i])
    node <- edge[row,1]
    edge_nr <- edge[row,2]
    row_node_daughter <- which(edge[,2]==node)

    print(row)
    print(node)
    print(edge_nr)
    print(row_node_daughter)

    edge[row_node_daughter,2] <- edge_nr
    edge <- edge[-row,]
    print(edge.length[-row])
    edge.length <- edge.length[-row]
    i <- i+1
  }

  edge_saved <- edge
  #Correcting for the fact that there are now "missing node numbers".
  n <- unique(edge[,1])
  nm <- matrix(0,length(n),3)
  nm[,1] <- n
  few_m <- few[few<max(edge[,1])]
  i <- 1
  for(i in 1:length(few_m)){

    pos <- which(nm[,1]>few_m[i])
    nm[pos,2] <- nm[pos,2]+1

    i <- i+1
  }
  nm[,3] <- nm[,1]-nm[,2]

  i <- 1
  for(i in 1:nrow(edge)){
    pos <- which(edge[i,1]==nm[,1])
    if(length(pos)!=0){
      edge[i,1] <- nm[pos,3]
    }
    pos <- which(edge[i,2]==nm[,1])
    if(length(pos)!=0){
      edge[i,2] <- nm[pos,3]
    }

    i <- i+1
  }

  #Nnode
  Nnode <- length(unique(edge[,1]))

  #Making the phylo object.
  x <- list()
  x$tip.label <- tip.label
  x$edge <- edge
  class(x$edge) <- "double"
  x$edge.length <- edge.length
  x$Nnode <- Nnode
  class(x) <- "phylo"

  y <- cbind(as.matrix(edge),as.matrix(edge.length))
  y <- y[order(y[,1]),]
  edge <- y[,1:2]
  edge.length <- t(y[,3])

  MMP_tree <- x

  #Saving the phylo object.
  print("Saving the phylo object")
  if(missing(tumorname)==TRUE){
    tumorname <- "No_name_declared_"
    print("You did not declare a name for the output file. Giving the name No_name_declared_")
  }
  write.xlsx(MMP_tree$edge,paste(tumorname,"_MMP_phylo.xlsx",sep=""),sheetName="edge")
  write.xlsx(MMP_tree$tip.label,paste(tumorname,"_MMP_phylo.xlsx",sep=""),sheetName="tip.label",append = TRUE,)
  write.xlsx(MMP_tree$Nnode,paste(tumorname,"_MMP_phylo.xlsx",sep=""),sheetName="Nnode",append = TRUE,)
  write.xlsx(MMP_tree$edge.length,paste(tumorname,"_MMP_phylo.xlsx",sep=""),sheetName="edge.length",append = TRUE,)

  print("Done")
  return(MMP_tree)
}
