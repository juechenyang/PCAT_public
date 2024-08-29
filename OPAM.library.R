OPAM.Projection.2 <- function(
    data.array,
    gene.names,
    n.cols,
    n.rows,
    weight = 0,
    statistic = "Kolmogorov-Smirnov",  # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
                                       # "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
                                       # "area.under.RES", or "Wilcoxon"
    gene.set,
    nperm = 200,
    correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{

    ES.vector <- vector(length=n.cols)
    NES.vector <- vector(length=n.cols)
    p.val.vector <- vector(length=n.cols)
    correl.vector <- vector(length=n.rows, mode="numeric")

# Compute ES score for signatures in each sample

#   print("Computing GSEA.....")
   phi <- array(0, c(n.cols, nperm))
   for (sample.index in 1:n.cols) {
      gene.list <- order(data.array[, sample.index], decreasing=T)            

      #      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" ")) 
      gene.set2 <- match(gene.set, gene.names)

      if (weight == 0) {
         correl.vector <- rep(1, n.rows)
      } else if (weight > 0) {
         if (correl.type == "rank") {
            correl.vector <- data.array[gene.list, sample.index]
         } else if (correl.type == "symm.rank") {
            correl.vector <- data.array[gene.list, sample.index]
            correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
                                    correl.vector,
                                    correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
         } else if (correl.type == "z.score") {
            x <- data.array[gene.list, sample.index]
            correl.vector <- (x - mean(x))/sd(x)
         }
      }
      GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
                             statistic = statistic, alpha = weight, correl.vector = correl.vector)
      ES.vector[sample.index] <- GSEA.results$ES

      if (nperm == 0) {
         NES.vector[sample.index] <- ES.vector[sample.index]
         p.val.vector[sample.index] <- 1
       } else {
         for (r in 1:nperm) {
            reshuffled.gene.labels <- sample(1:n.rows)
            if (weight == 0) {
               correl.vector <- rep(1, n.rows)
            } else if (weight > 0) {
                correl.vector <- data.array[reshuffled.gene.labels, sample.index]
            } 
            GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
                             statistic = statistic, alpha = weight, correl.vector = correl.vector)
            phi[sample.index, r] <- GSEA.results$ES
         }
         if (ES.vector[sample.index] >= 0) {
            pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
            if (length(pos.phi) == 0) pos.phi <- 0.5
            pos.m <- mean(pos.phi)
            NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
            s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
            p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
         } else {
            neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
            if (length(neg.phi) == 0) neg.phi <- 0.5 
            neg.m <- mean(neg.phi)
            NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
            s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
            p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
          }
       }
    }
    return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))

} # end of OPAM.Projection.2

OPAM.project.dataset.2 <- function( 
   input.ds,
   output.ds,
   gene.set.databases,
   gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
   sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
   weight              = 0.25,
   statistic           = "area.under.RES",
   output.score.type   = "ES",  # "ES" or "NES"
   nperm               = 200,  # number of random permutations for NES case
   combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
                                         # a single score. "combine.replace" combine *_UP and 
                                         # *_DN versions in a single score that replaces the individual
                                         # *_UP and *_DN versions. "combine.add" combine *_UP and 
                                         # *_DN versions in a single score and add it but keeping 
                                         # the individual *_UP and *_DN versions.
    correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------

   # Load libraries
   library(gtools)
   library(verification)
   library(RColorBrewer)

   # Read input dataset

   dataset <- MSIG.Gct2Frame(ds = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   gene.names <- dataset$row.names
   gene.descs <- dataset$descs
   sample.names <- dataset$names
   Ns <- length(m[1,])
   Ng <- length(m[,1])
   # temp <- strsplit(input.ds, split="/") # Extract input file name
   # s <- length(temp[[1]])
   # input.file.name <- temp[[1]][s]
   # temp <- strsplit(input.file.name, split=".gct")
   # input.file.prefix <-  temp[[1]][1]

    # Sample normalization

   if (sample.norm.type == "rank") {
      for (j in 1:Ns) {  # column rank normalization 
         m[,j] <- rank(m[,j], ties.method = "average")
      }
      m <- 10000*m/Ng
   } else if (sample.norm.type == "log.rank") {
      for (j in 1:Ns) {  # column rank normalization 
         m[,j] <- rank(m[,j], ties.method = "average")
      }
      m <- log(10000*m/Ng + exp(1))
   } else if (sample.norm.type == "log") {
      m[m < 1] <- 1
      m <- log(m + exp(1))
   }
  
   # Read gene set databases

   max.G <- 0
   max.N <- 0
   for (gsdb in gene.set.databases) {
      GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
      max.G <- max(max.G, max(GSDB$size.G))
      max.N <- max.N +  GSDB$N.gs
   }
   N.gs <- 0
   gs <- matrix("null", nrow=max.N, ncol=max.G)
   gs.names <- vector(length=max.N, mode="character")
   gs.descs <- vector(length=max.N, mode="character")
   size.G <- vector(length=max.N, mode="numeric")
   start <- 1
   for (gsdb in gene.set.databases) {
      GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
      N.gs <- GSDB$N.gs 
      gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
      gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
      size.G[start:(start + N.gs - 1)] <- GSDB$size.G
      gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
      start <- start + N.gs
   }
   N.gs <- max.N

   # Select desired gene sets

   if (gene.set.selection[1] != "ALL") {
      locs <- match(gene.set.selection, gs.names)
      N.gs <- sum(!is.na(locs))
      gs <- gs[locs,]
      gs.names <- gs.names[locs]
      gs.descs <- gs.descs[locs]
      size.G <- size.G[locs]
   }

   # Loop over gene sets

   score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
   for (gs.i in 1:N.gs) {
      gene.set <- gs[gs.i, 1:size.G[gs.i]]
      gene.overlap <- intersect(gene.set, gene.names)
      print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
      if (length(gene.overlap) == 0) { 
         score.matrix[gs.i, ] <- runif(Ns, min=1E-06, max=1.1E-06)
         next
      } else {
         gene.set.locs <- match(gene.overlap, gene.set)
         gene.names.locs <- match(gene.overlap, gene.names)
         msig <- m[gene.names.locs,]
         msig.names <- gene.names[gene.names.locs]
         if (output.score.type == "ES") {
            OPAM <- OPAM.Projection.2(data.array = m, gene.names = gene.names, n.cols = Ns, 
                                 n.rows = Ng, weight = weight, statistic = statistic,
                                 gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
            score.matrix[gs.i,] <- OPAM$ES.vector
         } else if (output.score.type == "NES") {
            OPAM <- OPAM.Projection.2(data.array = m, gene.names = gene.names, n.cols = Ns, 
                                 n.rows = Ng, weight = weight, statistic = statistic,
                                 gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
            score.matrix[gs.i,] <- OPAM$NES.vector
         }
      }
   }

   initial.up.entries <- 0
   final.up.entries <- 0
   initial.dn.entries <- 0
   final.dn.entries <- 0
   combined.entries <- 0
   other.entries <- 0

   if (combine.mode == "combine.off") {
      score.matrix.2 <- score.matrix
      gs.names.2 <- gs.names
      gs.descs.2 <- gs.descs
   } else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
      score.matrix.2 <- NULL
      gs.names.2 <- NULL
      gs.descs.2 <- NULL
      k <- 1
      for (i in 1:N.gs) {
            temp <- strsplit(gs.names[i], split="_") 
            body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
            suffix <- tail(temp[[1]], 1)
            print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
            if (suffix == "UP") {  # This is an "UP" gene set
               initial.up.entries <- initial.up.entries + 1
               target <- paste(body, "DN", sep="_")
               loc <- match(target, gs.names)            
               if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
                  score <- score.matrix[i,] - score.matrix[loc,]
                  score.matrix.2 <- rbind(score.matrix.2, score)
                  gs.names.2 <- c(gs.names.2, body)
                  gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
                  combined.entries <- combined.entries + 1
                  if (combine.mode == "combine.add") {  # also add the "UP entry
                     score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                     gs.names.2 <- c(gs.names.2, gs.names[i])
                     gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                     final.up.entries <- final.up.entries + 1
                  }
               } else { # did not find corresponding "DN" gene set: create "UP" entry
                  score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                  gs.names.2 <- c(gs.names.2, gs.names[i])
                  gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                  final.up.entries <- final.up.entries + 1
               }
            } else if (suffix == "DN") { # This is a "DN" gene set
               initial.dn.entries <- initial.dn.entries + 1
               target <- paste(body, "UP", sep="_")
               loc <- match(target, gs.names)            
               if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
                  score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                  gs.names.2 <- c(gs.names.2, gs.names[i])
                  gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                  final.dn.entries <- final.dn.entries + 1
               } else { # it found corresponding "UP" gene set
                  if (combine.mode == "combine.add") { # create "DN" entry
                     score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                     gs.names.2 <- c(gs.names.2, gs.names[i])
                     gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                     final.dn.entries <- final.dn.entries + 1
                  }
               }
            } else { # This is neither "UP nor "DN" gene set: create individual entry
               score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
               gs.names.2 <- c(gs.names.2, gs.names[i])
               gs.descs.2 <- c(gs.descs.2, gs.descs[i])
               other.entries <- other.entries + 1
            }
      } # end for loop over gene sets
      print(paste("initial.up.entries:", initial.up.entries))
      print(paste("final.up.entries:", final.up.entries))
      print(paste("initial.dn.entries:", initial.dn.entries))
      print(paste("final.dn.entries:", final.dn.entries))
      print(paste("other.entries:", other.entries))
      print(paste("combined.entries:", combined.entries))
      print(paste("total entries:", length(score.matrix.2[,1])))
   }            

   V.GCT <- data.frame(score.matrix.2)
   names(V.GCT) <- sample.names
   row.names(V.GCT) <- gs.names.2
   write.gct(gct.data.frame = V.GCT, descs = gs.descs.2, filename = output.ds)  

} # end of OPAM.project.dataset.2


OPAM.match.projection.to.phenotypes <-  function(
    input.ds,
    input.cls,
    results.dir,
    normalize.score = T,
    normalization.type = "zero.one",  # "zero.one", "z.score" or "r.z.score"
    markers.num=5,
    user.colors = NA,
    markers.metric = "ROC",   # "ROC" or "T.TEST"
    markers.file = NULL,
    sort.phenotypes = T,
    sort.decreasing = T,    # T = decreasing, F = increasing
    sort.expression = T,
    sort.decreasing.genes = T,
    legend = T,
    char.res = 1,
    only.up = F,
    cmap.type = 3,
    show.desc = T,
    row.norm = T)
  {

   library(gtools)
   library(verification)
   library(ROCR)
   library(MASS)
   library(RColorBrewer)
   library(heatmap.plus)

   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(m[1,])

   for (i in 1:length(m[,1])) {
      if (sd(m[i,]) == 0) {
         val <- m[i, 1]
         m[i,] <- m[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
      }
   }
   dim(m)
   sample.names <- dataset$names

   n.models <- length(m[,1])
   temp <- strsplit(input.ds, split="/") # Extract test file name
   s <- length(temp[[1]])
   test.file.name <- temp[[1]][s]
   temp <- strsplit(test.file.name, split=".gct")
   test.file.prefix <-  temp[[1]][1]
#   char.res <-  0.013 * n.models + 0.65

   # normalize scores

   if (normalize.score == T) {
     if (normalization.type == "zero.one") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
         }
     } else if (normalization.type == "z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
         }
     } else if (normalization.type == "r.z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
          }
     }         
   }

   CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   if (is.vector(cls.labels)) {
      n.phen <- 1
   } else {
      n.phen <- length(cls.labels[,1])
   }
   if (!is.na(user.colors)) {
      c.test <- user.colors
    } else {
      if (!is.null(CLS$col.phen)) {
         c.test <- CLS$col.phen
      } else {
         c.test <- c(brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"))
      }
    }


   if (!is.null(CLS$phen.names)) {
      phen.names <- CLS$phen.names
#      if (is.vector(cls.list)) {
#         cls.phen <- paste(phen.names, cls.phen, collapse="_")
#      } else {
#         for (i in 1:length(cls.phen)) {
#            for (j in 1:length(cls.phen[[i]])) {
#               cls.phen[[i]][j] <- paste(phen.names[i], cls.phen[[i]][j], collapse="_")
#            }
#         }
#      }
   } else {
      phen.names <- "NA"
   }

   cls.phen.index <- unlist(cls.phen)
   cls.phen.colors <- c.test[1:length(cls.phen.index)]

   n.classes <- vector(length=n.phen, mode="numeric")
   if (n.phen == 1) {
      max.classes <- length(cls.phen)
      n.classes[1] <- max.classes
   } else {
     max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
     for (i in 1:n.phen) {
       n.classes[i] <- length(cls.phen[[i]])
     }
   }

   x <- rbind(sample.names, cls.list, cls.labels)
   print("before loop")
   print(x)
   print(cls.phen)
   print(phen.names)

   filename <- paste(results.dir, test.file.prefix, ".PHEN.MARKERS.", markers.metric, ".pdf", sep="")
   pdf(file=filename, height = 10, width = 10)

   # Loop over phenotypes

   for (k.phen in 1:n.phen) {


      if (is.vector(cls.labels)) {
         k.phen.labels <- cls.labels
         k.phen.list <- cls.list
      } else {
         k.phen.labels <- as.vector(cls.labels[k.phen,])
         k.phen.list <- as.vector(cls.list[k.phen,])
      }

      # Sort according to current phenotype

      if(sort.expression == T) {
         phen.index <- order(k.phen.labels, decreasing=sort.decreasing)
      } else {
         phen.index <- seq(1, length(k.phen.labels))
      }
      if (is.vector(cls.labels)) {
         cls.labels2 <- cls.labels[phen.index]
         cls.list2 <- cls.list[phen.index]
      } else {
         cls.labels2 <- cls.labels[, phen.index]
         cls.list2 <- cls.list[, phen.index]
      }
      k.phen.labels <- k.phen.labels[phen.index]
      k.phen.list <- k.phen.list[phen.index]
      sample.names2 <- sample.names[phen.index]
      m2 <- m[, phen.index]

   x <- rbind(sample.names2, cls.list2, cls.labels2)
   print(paste("inside loop phen=", k.phen))
   print(x)
   print(cls.phen)
   print(phen.names)

      # Markers for each class

      if (is.vector(cls.labels2)) {
         classes <- unique(cls.list2)
      } else {
         classes <- unique(cls.list2[k.phen, ])
      }
      if (length(classes) > 2) {
         k.only.up <- T
      } else {
         k.only.up <- only.up
      }

      if(length(classes) == 2) classes <- classes[1]
      markers <- NULL
      markers.descs <- NULL
      metric.list <- NULL
      p.val.list <- NULL
      k.class <- NULL
      for (k in classes) {
         if (is.vector(cls.labels2)) {
            bin.class <- ifelse(cls.list2 == k, 0, 1)
         } else {
            bin.class <- ifelse(cls.list2[k.phen, ] == k, 0, 1)
         }
         if (markers.metric == "T.TEST") {
            metric <- vector(length=n.models, mode="numeric")
            p.val <- vector(length=n.models, mode="numeric")
            for (i in 1:n.models) {
               temp <- split(m2[i, ], bin.class)
               x <- temp[[1]]
               y <- temp[[2]]
               metric[i] <- signif(t.test(x=x, y=y)$statistic, digits=3)
               p.val[i] <- signif(t.test(x=x, y=y)$p.value, digits=3)
            }
         } else if (markers.metric == "ROC") {
            bin.class <- ifelse(bin.class == 1, 0, 1)
            metric <- vector(length=n.models, mode="numeric")
            p.val <- vector(length=n.models, mode="numeric")
            for (i in 1:n.models) {
               m.score <- m2[i,]
               m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
               if (length(table(bin.class)) > 1) {
                  perf.auc <- roc.area(bin.class, m.score.norm)
                  metric[i] <- signif(perf.auc$A, digits=3)
                  p.val[i] <- signif(perf.auc$p.value, digits=3)
               } else {
                  metric[i] <- 1
                  p.val[i] <- 1
               }
            }
         } else if (markers.metric == "MEAN.DIFF") {
            bin.class <- ifelse(bin.class == 1, 0, 1)
            metric <- vector(length=n.models, mode="numeric")
            p.val <- vector(length=n.models, mode="numeric")
            for (i in 1:n.models) {
               temp <- split(m2[i, ], bin.class)
               x <- temp[[1]]
               y <- temp[[2]]
               metric[i] <- signif(mean(x) - mean(y), digits=3)
               p.val[i] <- signif(t.test(x=x, y=y)$p.value, digits=3)
            }
         }

         if (is.na(sort.decreasing.genes)) {
            metric.order <- seq(1, length(metric))
         } else {
            metric.order <- order(metric, decreasing=sort.decreasing.genes)
         }
         if (only.up == TRUE) {
            k.markers.num <- ifelse(markers.num > n.models, n.models, markers.num)

#            if (length(classes) == 2) {
#               k.markers.num <- ifelse(markers.num > n.models, n.models, markers.num)
#            } else {
#               k.markers.num <- ifelse(length(classes)*markers.num > n.models, 
#                                               floor(n.models/length(classes)), markers.num)
#            }
            markers <- c(markers, model.names[metric.order][1:k.markers.num])
            markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num])
            metric.list <- c(metric.list, metric[metric.order][1:k.markers.num])
            p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num])
            k.class <- c(k.class, rep(k, k.markers.num))
         } else {
            k.markers.num <- ifelse(length(classes)*markers.num > n.models, floor(n.models/length(classes)), 
                                                                          markers.num)
            markers <- c(markers, model.names[metric.order][1:k.markers.num],
                         model.names[metric.order][(length(model.names) - k.markers.num +1):length(model.names)])
            markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num],
                               model.descs[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
            metric.list <- c(metric.list, metric[metric.order][1:k.markers.num],
                             metric[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
            p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num],
                            p.val[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
            k.class <- c(k.class, rep(k, k.markers.num), rep(paste("not", k), k.markers.num))
         }
      }

      V3 <- m2[markers,]
      print(V3)
      print(markers)

      if (show.desc == T) {
         model.descs2 <- paste(metric.list, p.val.list, k.class, markers.descs)
      } else {
         model.descs2 <- paste(metric.list, p.val.list)
      }
      height <- ifelse(length(markers) + n.phen >= 9, 10, (length(markers) + n.phen)*0.44 + 5)
#      char.res <-  0.0085 * length(markers) + 0.65


      # Sort markers inside each phenotype class

      if(sort.expression == T) {
         for (j in unique(k.phen.labels)) {
            V4 <- V3[ , k.phen.labels == j]
            sn <- sample.names2[k.phen.labels == j]
            if (is.vector(cls.labels)) {
               clab <- cls.labels2[k.phen.labels == j]
               clis <- cls.list2[k.phen.labels == j]
            } else {
               clab <- cls.labels2[, k.phen.labels == j]
               clis <- cls.list2[, k.phen.labels == j]
            }
            l.phen <- sum(k.phen.labels == j)
            if (l.phen > 1) {
               dist.matrix <- dist(t(V4))
               HC <- hclust(dist.matrix, method="complete")
               HC.order <- HC$order
               V4 <- V4[ , HC.order]
               sn <- sn[HC.order]
               if (is.vector(cls.labels2)) {
                  clab <- clab[HC.order]
                  clis <- clis[HC.order]
               } else {
                  clab <- clab[, HC.order]
                  clis <- clis[, HC.order]
               }
            }
            V3[ , k.phen.labels == j] <- V4
            sample.names2[k.phen.labels == j] <- sn
            if (is.vector(cls.labels2)) {
               cls.labels2[k.phen.labels == j] <- clab
               cls.list2[k.phen.labels == j] <- clis
            } else {
               cls.labels2[, k.phen.labels == j] <- clab
               cls.list2[, k.phen.labels == j] <- clis
            }
         }
   }
   x <- rbind(sample.names2, cls.list2, cls.labels2)
   print(paste("inside loop after in-class sort phen=", k.phen))
   print(x)
   print(cls.phen)
   print(phen.names)

     # Recompute cls.phen and cls.labels2 as order may have changed

     cls.phen2 <- NULL
     if (is.vector(cls.labels2)) {
        classes <- unique(cls.list2)
        cls.phen2 <- classes
        cls.labels2 <- match(cls.list2, cls.phen2)
      } else {
         for (kk in 1:length(cls.list2[, 1])) {
            classes <- unique(cls.list2[kk,])
            cls.phen2[[kk]] <- classes
            cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
         }
      }

   x <- rbind(sample.names2, cls.list2, cls.labels2)
   print(paste("inside loop after cls.phen renorm phen=", k.phen))
   print(cls.phen2)
   print(phen.names)


     library(gmodels)
     if (!is.vector(cls.labels2)) {
        if (sort.phenotypes == T) {
           phen.score <- vector(length=n.phen, mode="numeric")
           for (k.lab in 1:n.phen) {
              tab <- table(as.vector(cls.list2[k.lab,]), k.phen.list)
              print(tab)
#              phen.score[k.lab] <- 1 - chisq.test(tab)$p.value
#              phen.score[k.lab] <- 1 - fisher.test(tab)$p.value
              if ((length(tab[,1]) > 1) && (length(tab[1,]) > 1)) { 
                 CT <- CrossTable(tab, chisq=T)
                 phen.score[k.lab] <- CT$chisq$p.value
                 print(phen.score[k.lab])
              } else {
                 phen.score[k.lab] <- 0.50
                 print(phen.score[k.lab])
              }
           }
           phen.order <- order(phen.score, decreasing= T)
           print(phen.order)
           cls.labels2 <- cls.labels2[phen.order,]
           cls.phen2 <- cls.phen2[phen.order]
           phen.names2 <- phen.names[phen.order]
           main.string <- paste(test.file.prefix, " - ", phen.names2[n.phen], markers.metric, " order")
        } else {
           phen.names2 <- phen.names
           main.string <- paste(test.file.prefix, " - ", phen.names2[k.phen], markers.metric, " order")
        }
     } else {
        phen.names2 <- phen.names[1]
        main.string <- paste(test.file.prefix, " - ", phen.names2, markers.metric, " order")
     }

#     windows(width=15, height=height)


   x <- rbind(sample.names2, cls.list2, cls.labels2)
   print(paste("inside loop after phen sort before figure phen=", k.phen))
   print(x)
   print(cls.phen2)
   print(phen.names2)

   phen.list <- unlist(cls.phen2)
   colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
   
   print(rbind(phen.list, colors.list))

   if (show.desc == T) {
      markers <- paste(markers, seq(1, length(markers)), sep="_")
   }

   MSIG.HeatMapPlot.7(V = V3, row.names = markers,
                      row.names2 = model.descs2, col.labels = cls.labels2, 
                      col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names2,
                      col.names = sample.names2, main = main.string, xlab="  ", ylab="  ", 
                      row.norm = row.norm,  
                      cmap.type = cmap.type, char.rescale = char.res,  legend=legend)

     V3 <- data.frame(V3)
     colnames(V3) <- sample.names2
     row.names(V3) <- markers

     if (!is.null(markers.file)) {
        write.gct(gct.data.frame = V3, descs = model.descs2, filename = markers.file)  
     }

   } # end loop over phenotypes

   dev.off()
    
}

OPAM.sort.projection.by.score.2 <- function(
    input.ds,
    input.cls,
    results.dir,
    normalize.score = T,
    normalization.type = "zero.one",
    model,
    target.phen = NA,
    target.class = NA,
    user.colors = NA,
    decreasing.order = T,
    output.dataset = NA,
    char.rescale = 1,
    cmap.type = 3,
    row.norm = T)
  {

   library(gtools)
   library(verification)
   library(ROCR)
   library(MASS)
   library(RColorBrewer)
   library(heatmap.plus)

   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(m[1,])
   dim(m)
   sample.names <- dataset$names

   n.models <- length(m[,1])
   temp <- strsplit(input.ds, split="/") # Extract test file name
   s <- length(temp[[1]])
   test.file.name <- temp[[1]][s]
   temp <- strsplit(test.file.name, split=".gct")
   test.file.prefix <-  temp[[1]][1]
   char.res <-  0.013 * n.models + 0.65

   # normalize scores

   if (normalize.score == T) {
     if (normalization.type == "zero.one") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
         }
     } else if (normalization.type == "z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
         }
     } else if (normalization.type == "r.z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
          }
     }         
   }

   CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   if (is.vector(cls.labels)) {
      n.phen <- 1
   } else {
      n.phen <- length(cls.labels[,1])
   }
   if (!is.na(user.colors)) {
      c.test <- user.colors
    } else {
      if (!is.null(CLS$col.phen)) {
         c.test <- CLS$col.phen
      } else {
         c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"))
      }
    }


   if (!is.null(CLS$phen.names)) {
      phen.names <- CLS$phen.names
   } else {
      phen.names <- "NA"
   }

   cls.phen.index <- unlist(cls.phen)
   cls.phen.colors <- c.test[1:length(cls.phen.index)]

   n.classes <- vector(length=n.phen, mode="numeric")
   if (n.phen == 1) {
      max.classes <- length(cls.phen)
      n.classes[1] <- max.classes
   } else {
     max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
     for (i in 1:n.phen) {
       n.classes[i] <- length(cls.phen[[i]])
     }
   }

   filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
   pdf(file=paste(filename, ".pdf", sep=""), height = 8.5, width = 11)
#   windows(width=12, height=8)

   loc <- match(model, model.names)
   print(c("loc:", loc))
   s.order <- order(m[loc,], decreasing = decreasing.order)
   m2 <- m[, s.order]

   sample.names2 <- sample.names[s.order]

   if (is.vector(cls.labels)) {
      cls.labels2 <- cls.labels[s.order]
      cls.list2 <- cls.list[s.order]
   } else {
      cls.labels2 <- cls.labels[, s.order]
      cls.list2 <- cls.list[, s.order]
   }
      # Recompute cls.phen and cls.labels2 as order may have changed

     cls.phen2 <- NULL
     if (is.vector(cls.labels)) {
        classes <- unique(cls.list2)
        cls.phen2 <- classes
        cls.labels2 <- match(cls.list2, cls.phen2)
      } else {
         for (kk in 1:length(cls.list2[, 1])) {
            classes <- unique(cls.list2[kk,])
            cls.phen2[[kk]] <- classes
            cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
         }
   }


   correl <- cor(t(m2))[, loc]
   m.order <- order(correl, decreasing=T)
   correl2 <- correl[m.order]
   m2 <- m2[m.order,]
   model.names2 <- model.names[m.order]
   model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
   phen.list <- unlist(cls.phen2)
   colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
 
   if (!is.na(target.phen)) {
       bin.class <- ifelse(cls.list2[target.phen,] == target.class, 1, 0)
   } else {
       bin.class <- ifelse(cls.list2[1,] == cls.list2[1,1], 1, 0)
   }
   for (i in 1:n.models) {
      m.score <- m2[i,]
      m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
      perf.auc <- roc.area(bin.class, m.score.norm)
      print(paste("ROC=", signif(perf.auc$A, digits=3), " p-val=", signif(perf.auc$p.value, digits=3))) 
      model.descs2[i] <- paste(signif(perf.auc$A, digits=3), " (", signif(perf.auc$p.value, digits=3), ")")
   }
      
   MSIG.HeatMapPlot.7(V = m2, row.names = model.names2,
                      row.names2 = model.descs2, col.labels = cls.labels2, 
                      col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
                      col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = row.norm,  
                      cmap.type = cmap.type, char.rescale = char.rescale,  legend=T)

   dev.off()

   if (!is.na(output.dataset)) {
      V.GCT <- m2
      colnames(V.GCT) <- sample.names2
      row.names(V.GCT) <- model.names2
      write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename =output.dataset)  
   }
    
 }



MSIG.HeatMapPlot.7 <- function(
V, 
row.names = "NA",
row.names2 = "NA", 
col.labels = "NA",
col.labels2 = "NA", 
col.classes = "NA", 
phen.cmap = "NA", 
col.names = "NA",
phen.names = "NA",                               
main = " ", 
sub = " ", 
xlab=" ", 
ylab=" ",
row.norm = TRUE,
char.rescale = 0.85,                               
cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
max.v = "NA",
legend = T)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

       n.rows <- length(V[,1])
       n.cols <- length(V[1,])
       V1 <- matrix(0, nrow=n.rows, ncol=n.cols)

#       if ((cmap.type == 5) | (cmap.type == 3)) {
              if (cmap.type == 5) {
          row.norm <- F
       }

       if (row.norm == TRUE) {
          row.mean <- apply(V, MARGIN=1, FUN=mean)
          row.sd <- apply(V, MARGIN=1, FUN=sd)
          row.n <- length(V[,1])
          for (i in 1:n.rows) {
	     if (row.sd[i] == 0) {
    	         V1[i,] <- 0
             } else {
	         V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
             }
             V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
             V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
          }
        } else {
          V1 <- V
        }

        if (cmap.type == 1) { 
             mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
                        "#FF9DB0", "#FF7080", 
                        "#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
                                                         # pre-gene cluster, original pinkogram color map
        } else if (cmap.type == 2) {
           violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
           mycol <- rev(violet.palette(20))

#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
        } else if (cmap.type == 6) {
             mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
                        "#7A0177", "#49006A")
        } else if (cmap.type == 7) {
             mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
                        "#A63603", "#7F2704")
        } else if (cmap.type == 8) {
            mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
                       "#006D2C", "#00441B")
        } else if (cmap.type == 9) {
            mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
                       "#08519C", "#08306B")
        } else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
            mycol <- vector(length=512, mode = "numeric")

            for (k in 1:256) {
               mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
            }
            for (k in 257:512) {
               mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
            }
            mycol <- rev(mycol)
          }

       ncolors <- length(mycol)

       if (cmap.type == 5) {
           if (max.v == "NA") {
              max.v <- max(max(V1), -min(V1))
            }
           V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))

       } else {
           V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
        }

        if (col.labels[1] == "NA") {      
           heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
           heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
           tot.cols <- ncolors
           if (legend == T) {
              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
           } else {
              nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
           }
           par(mar = c(3, 16, 3, 16))
           mycol <- c(mycol, phen.cmap[1:length(col.classes)])
           image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
              main=main, sub = sub, xlab= xlab, ylab=ylab)
           n.rows.phen <- 0
         } else {
           tot.cols <- ncolors
           if (is.vector(col.labels)) {
              heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
              heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
              n.rows.phen <- 1
              heatm[n.rows + 1,] <- tot.cols + col.labels
              cols.row <- length(unique(col.labels))
              tot.cols <- tot.cols + cols.row
              phen.cmap <- phen.cmap[1:cols.row]
            } else {
              n.rows.phen <- length(col.labels[,1])
              cols.row <- vector(length=n.rows.phen, mode = "numeric")
              heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
              heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
              for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
                 heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
                 cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
                 tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))

               }
              phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
            }
           if (legend == T) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
              nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
           } else {
              nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
           }
           par(mar = c(3, 16, 3, 16))
           mycol <- c(mycol, phen.cmap)
           image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
                 main=main, sub = sub, xlab= xlab, ylab=ylab)
         }

# Add lines to separate phenotypes or subgroups

       if (col.labels2[1] != "NA") {
          groups <-  split(col.labels2, col.labels2)
          len.vec <- lapply(groups, length)
          plot.div <- c(0.51, cumsum(len.vec) + 0.5)
          for (i in plot.div) {
             lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
          }
          lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
          lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
                cex = 0.9, col = "black")
          lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
                cex = 0.9, col = "black")
        }
       if (row.names[1] != "NA") {
            numC <- nchar(row.names)
            size.row.char <- char.rescale*25/(n.rows + 20)
            for (i in 1:n.rows) {
               row.names[i] <- substr(row.names[i], 1, 40)
               row.names[i] <- paste(row.names[i], " ", sep="")
            }
            if (phen.names[1] == "NA") {
               head.names <- paste("Class", seq(n.rows.phen, 1, -1))
             } else {
               head.names <- as.character(rev(phen.names))
             }
            row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
            axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
                 font.axis=2, line=-1)
        }

       if (row.names2[1] != "NA") {
            numC <- nchar(row.names2)
            size.row.char <- char.rescale*25/(n.rows + 20)
            for (i in 1:n.rows) {
               row.names2[i] <- substr(row.names2[i], 1, 40)
               row.names2[i] <- paste(" ", row.names2[i], sep="")
            }
            row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "   ")
            axis(4, at=1:(n.rows + 1), labels=row.names2, adj= 0.5, tick=FALSE, las = 1, 
                 cex.axis=size.row.char, font.axis=2, line=-1)
        }

        if (col.names[1] != "NA") {
          size.col.char <- char.rescale*20/(n.cols + 25)
          axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
        }

      # Phenotype Legend 

#      print("--------------------------------------------------------------------------------------------")
       if (legend == T) {
          leg.txt <- NULL
          p.vec <- NULL
          c.vec <- NULL
          c2.vec <- NULL
          ind <- 1
          par(mar = c(0, 0, 0, 0))
          plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
          for (i in 1:n.rows.phen) {  
            if (is.vector(col.labels)) {
                phen.v <- as.character(col.classes)
            } else {
                phen.v <- as.character(col.classes[[i]])
            }
            p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
            leg.txt <- c(p.name, phen.v)  
            p.vec <-  rep(22, cols.row[i] + 1)
            c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
            c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
            ind <- ind + cols.row[i]
            offset <- 0.07
            legend(x=0, y= 1 - offset*i, 
              horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
              pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
          }
        }
       
       # Color map legend

#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
       
       par(mar = c(2, 12, 2, 12))
       num.v <- 20
          range.v <- range(V2)
          incr <-  (range.v[1] - range.v[2])/(num.v - 1)
          heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
          image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
                main=" ", sub = " ", xlab= ylab, ylab=xlab)
          range.v <- range(V1)
          incr <-  (range.v[1] - range.v[2])/(num.v - 1)
          heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
#          print(c("heatm.v2=", heatm.v2))
          axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
              
	return()

     }

MSIG.Gct2Frame <- function(ds = NULL) {
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   # ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
   descs <- ds[,1]
   ds <- ds[-1]
   row.names <- row.names(ds)
   names <- names(ds)
   return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

Read.GeneSets.db <- function(
   gs.db,
   thres.min = 2,
   thres.max = 2000,
   gene.names = NULL)
  {

   temp <- readLines(gs.db)
   max.Ng <- length(temp)
   temp.size.G <- vector(length = max.Ng, mode = "numeric") 
   for (i in 1:max.Ng) {
      temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
   }
   max.size.G <- max(temp.size.G)      
   gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
   temp.names <- vector(length = max.Ng, mode = "character")
   temp.desc <- vector(length = max.Ng, mode = "character")
   gs.count <- 1
   for (i in 1:max.Ng) {
      gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
      gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
      gene.set.name <- gs.line[1] 
      gene.set.desc <- gs.line[2] 
      gene.set.tags <- vector(length = gene.set.size, mode = "character")
      for (j in 1:gene.set.size) {
         gene.set.tags[j] <- gs.line[j + 2]
      }
      if (is.null(gene.names)) {
         existing.set <- rep(TRUE, length(gene.set.tags))
      } else {
         existing.set <- is.element(gene.set.tags, gene.names)
      }
      set.size <- length(existing.set[existing.set == T])
      if ((set.size < thres.min) || (set.size > thres.max)) next
      temp.size.G[gs.count] <- set.size
      gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
      temp.names[gs.count] <- gene.set.name
      temp.desc[gs.count] <- gene.set.desc
      gs.count <- gs.count + 1
    }
   Ng <- gs.count - 1
   gs.names <- vector(length = Ng, mode = "character")
   gs.desc <- vector(length = Ng, mode = "character")
   size.G <- vector(length = Ng, mode = "numeric") 
   
   gs.names <- temp.names[1:Ng]
   gs.desc <- temp.desc[1:Ng]
   size.G <- temp.size.G[1:Ng]
   
   return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
 }

OPAM.Projection <- function(
    data.array,
    gene.names,
    n.cols,
    n.rows,
    weight = 0,
    statistic = "Kolmogorov-Smirnov",  # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
                                       # "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
                                       # "area.under.RES", or "Wilcoxon"
    gene.set,
    nperm = 200) {

    ES.vector <- vector(length=n.cols)
    NES.vector <- vector(length=n.cols)
    p.val.vector <- vector(length=n.cols)
    correl.vector <- vector(length=n.rows, mode="numeric")

# Compute ES score for signatures in each sample

#   print("Computing GSEA.....")
   phi <- array(0, c(n.cols, nperm))
   for (sample.index in 1:n.cols) {
      gene.list <- order(data.array[, sample.index], decreasing=T)            

      #      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" ")) 
      gene.set2 <- match(gene.set, gene.names)

      if (weight == 0) {
         correl.vector <- rep(1, n.rows)
      } else if (weight > 0) {
         correl.vector <- data.array[gene.list, sample.index]
      }
      GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
                             statistic = statistic, alpha = weight, correl.vector = correl.vector)
      ES.vector[sample.index] <- GSEA.results$ES

      if (nperm == 0) {
         NES.vector[sample.index] <- ES.vector[sample.index]
         p.val.vector[sample.index] <- 1
       } else {
         for (r in 1:nperm) {
            reshuffled.gene.labels <- sample(1:n.rows)
            if (weight == 0) {
               correl.vector <- rep(1, n.rows)
            } else if (weight > 0) {
                correl.vector <- data.array[reshuffled.gene.labels, sample.index]
            } 
            GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
                             statistic = statistic, alpha = weight, correl.vector = correl.vector)
            phi[sample.index, r] <- GSEA.results$ES
         }
         if (ES.vector[sample.index] >= 0) {
            pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
            if (length(pos.phi) == 0) pos.phi <- 0.5
            pos.m <- mean(pos.phi)
            NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
            s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
            p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
         } else {
            neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
            if (length(neg.phi) == 0) neg.phi <- 0.5 
            neg.m <- mean(neg.phi)
            NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
            s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
            p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
          }
       }
    }
    return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))

} # end of OPAM.Projection

write.gct <- function(gct.data.frame, descs = "", filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")

    names <- names(gct.data.frame)
    cat("\t", names[1], file = f, append = TRUE, sep = "")

    if (length(names) > 1) {
       for (j in 2:length(names)) {
           cat("\t", names[j], file = f, append = TRUE, sep = "")
       }
     }
    cat("\n", file = f, append = TRUE, sep = "\t")

    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
    m[, 1] <- row.names(gct.data.frame)
    if (length(descs) > 1) {
        m[, 2] <- descs
    } else {
        m[, 2] <- row.names(gct.data.frame)
    }
    index <- 3
    for (i in 1:dim(gct.data.frame)[2]) {
        m[, index] <- gct.data.frame[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)

}
MSIG.ReadPhenFile <- function(file = "NULL") { 
#
# Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
#  (numeric and character) for the samples in a gene expression file (RES or GCT format)
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

      cls.cont <- readLines(file)
      num.lines <- length(cls.cont)
      temp <- unlist(strsplit(cls.cont[[1]], " "))
      if (length(temp) == 3) {
         phen.names <- NULL
         col.phen <- NULL
       } else {
         l.phen.names <- match("phen.names:", temp)
         l.col.phen <- match("col.phen:", temp)
         phen.names <- temp[(l.phen.names + 1):(l.col.phen - 1)]
         col.phen <- temp[(l.col.phen + 1):length(temp)]
       }
      temp <- unlist(strsplit(cls.cont[[2]], " "))
      phen.list <- temp[2:length(temp)]

      for (k in 1:(num.lines - 2)) {
        temp <- unlist(strsplit(cls.cont[[k + 2]], " "))
        if (k == 1) {
           len <- length(temp)
           class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
           class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
           phen <- NULL
        }
        class.list[k, ] <- temp
        classes <- unique(temp)
        class.v[k, ] <- match(temp, classes)
        phen[[k]] <- classes
      }
      if (num.lines == 3) {
         class.list <- as.vector(class.list)
         class.v <- as.vector(class.v)
         phen <- unlist(phen)
       }
      return(list(phen.list = phen.list, phen = phen, phen.names = phen.names, col.phen = col.phen,
                  class.v = class.v, class.list = class.list))
}



GSEA.EnrichmentScore5 <- function(
#
# Computes the weighted GSEA score of gene.set in gene.list.
#
# This version supports multiple statistics:
#  "Kolmogorov-Smirnov" "Cramer-von-Mises" "Anderson-Darling" "Zhang_A" "Zhang_C" "Zhang_K" and "area.under.RES"
#
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
# necessary to input the correlation vector with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: "Kolmogorov-Smirnov" "Cramer-von-Mises" "Anderson-Darling" "Zhang_A" "Zhang_C" "Zhang_K", "area.under.RES", or "Wilcoxon"
#  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (e.g. real number between -1 and +1 for KS) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2008 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
gene.list,              # The ordered gene list (e.g. integers indicating the original position in the input dataset)  
gene.set,               # A gene set (e.g. integers indicating the location of those genes in the input dataset) 
statistic = "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises", "Anderson-Darling", "Zhang_A", "Zhang_C",
                                  #  "Zhang_K", "area.under.RES", or "Wilcoxon"
alpha = 1,              # The weight exponent
correl.vector = NULL    # A correlation vector of genes with phenotype or other appropriate weight
) {  

   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 
   orig.correl.vector <- correl.vector
   if (alpha == 0) correl.vector <- rep(1, N)   # unweighted case
   correl.vector <- abs(correl.vector)^alpha
   sum.correl  <- sum(correl.vector[tag.indicator == 1])
   P0 <- no.tag.indicator / Nm
   F0 <- cumsum(P0)
   Pn <- tag.indicator * correl.vector / sum.correl
   Fn <- cumsum(Pn)
   if (statistic == "Kolmogorov-Smirnov") {
      RES <- Fn - F0
      max.ES <- max(RES)
      min.ES <- min(RES)
      if (max.ES > - min.ES) {
         ES <- signif(max.ES, digits = 5)
         arg.ES <- which.max(RES)
      } else {
         ES <- signif(min.ES, digits=5)
         arg.ES <- which.min(RES)
      }
      return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
   } else if (statistic == "Cramer-von-Mises") {
            # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
            # Modified to separate positive and negative enrichment parts
      RES <- Fn - F0
      X <- RES^2 * P0
      X_p <- X[RES >= 0]
      X_n <- X[RES < 0]
      ES_p <- sqrt(sum(X_p)/N)
      ES_n <- sqrt(sum(X_n)/N)
      if (ES_p > ES_n) {
         ES <- signif(ES_p, digits = 5)
         arg.ES <- which.min(abs(X - max(X_p)))
      } else {
         ES <- - signif(ES_n, digits=5)
         arg.ES <- which.min(abs(X - max(X_n)))
      }      
      return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
   } else if (statistic == "Anderson-Darling") {
           # Based on Choulakian et al Canadian J. of Statistics 22, 125 (1994)
           # Modified to separate positive and negative enrichment parts
      RES <- Fn - F0
      F0_factor <- ifelse(F0 < 1/Nm | F0 > (Nm - 1)/Nm, rep(1, N), F0 * (1 - F0))
      X <- RES^2 * P0 / F0_factor
      X_p <- X[RES >= 0]
      X_n <- X[RES < 0]
      ES_p <- sqrt(sum(X_p)/N)
      ES_n <- sqrt(sum(X_n)/N)
      if (ES_p > ES_n) {
         ES <- signif(ES_p, digits = 5)
         arg.ES <- which.min(abs(X - max(X_p)))
      } else {
         ES <- - signif(ES_n, digits=5)
         arg.ES <- which.min(abs(X - max(X_n)))
      }      
      return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
   } else if (statistic == "Zhang_A") {
           # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
           # Modified to separate positive and negative enrichment parts
      RES <- Fn - F0
      Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn * log(Fn/F0))
      Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl - 1)/sum.correl, 0, (1 - Fn) * log( (1 - Fn) / (1 - F0) ))
      Fn_factor <- ifelse(Fn < 1/sum.correl | Fn > (sum.correl - 1)/sum.correl, rep(1, N), Fn * (1 - Fn))
      G <- (Fact1 + Fact2) * Pn / Fn_factor 
      G_p <- G[RES >= 0]
      G_n <- G[RES < 0]
      ES_p <- sum(G_p)/N
      ES_n <- sum(G_n)/N
      if (ES_p > ES_n) {
         ES <- signif(ES_p, digits = 5)
         arg.ES <- which.min(abs(G - max(G_p)))
      } else {
         ES <- - signif(ES_n, digits=5)
         arg.ES <- which.min(abs(G - max(G_n)))
      }      
      return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
   } else if (statistic == "Zhang_C") {
           # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
           # Modified to separate positive and negative enrichment parts
      RES <- Fn - F0
      Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn * log(Fn/F0))
      Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl - 1)/sum.correl, 0, (1 - Fn) * log( (1 - Fn) / (1 - F0) ))
      F0_factor <- ifelse(F0 < 1/Nm | F0 > (Nm - 1)/Nm, rep(1, N), F0 * (1 - F0))
      G <- (Fact1 + Fact2) * P0 / F0_factor
      G_p <- G[RES >= 0]
      G_n <- G[RES < 0]
      ES_p <- sum(G_p)/N
      ES_n <- sum(G_n)/N
      if (ES_p > ES_n) {
         ES <- signif(ES_p, digits = 5)
         arg.ES <- which.min(abs(G - max(G_p)))
       } else {
         ES <- - signif(ES_n, digits=5)
         arg.ES <- which.min(abs(G - max(G_n)))
      }      
      return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
   } else if (statistic == "Zhang_K") {
           # Based on Zhang, J. R. Statist. Soc. B 64, Part2, 281 (2002)
           # Modified to separate positive and negative enrichment parts
      RES <- Fn - F0
      Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn * log(Fn/F0))
      Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl - 1)/sum.correl, 0, (1 - Fn) * log( (1 - Fn) / (1 - F0) ))
      G <- Fact1 + Fact2
      G_p <- G[RES >= 0]
      G_n <- G[RES < 0]
      ES_p <- max(G_p)
      ES_n <- max(G_n)
      if (ES_p > ES_n) {
         ES <- signif(ES_p, digits = 5)
         arg.ES <- which.min(abs(G - ES_p))
      } else {
         ES <- - signif(ES_n, digits=5)
         arg.ES <- which.min(abs(G - ES_n))
      }
      return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
   } else if (statistic == "area.under.RES") {
     # Area under the RES score
      RES <- Fn - F0
      max.ES <- max(RES)
      min.ES <- min(RES)
      if (max.ES > - min.ES) {
         arg.ES <- which.max(RES)
      } else {
         arg.ES <- which.min(RES)
      }
      ES <- sum(RES)
      return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
   } else if (statistic == "Wilcoxon") {
     # Wilcoxon test score
      library(exactRankTests)
      seq.index <- seq(1, N)
      gene.set.ranks <- seq.index[tag.indicator == 1]
      gene.set.comp.ranks <- seq.index[tag.indicator == 0]
      W <- wilcox.exact(x=gene.set.ranks, y =gene.set.comp.ranks, alternative = "two.sided", mu = 0,
                        paired = FALSE, exact = F, conf.int = T, conf.level = 0.95)
      ES <- log(1/W$p.value)
      return(list(ES = ES, arg.ES = NULL, RES = NULL, indicator = tag.indicator))
    }
 }

OPAM.sort.projection.by.score <- function(
    input.ds,
    input.cls,
    results.dir,
    normalize.score = T,
    normalization.type = "zero.one",
    model,
    user.colors = NA,
    decreasing.order = T)
  {

   library(gtools)
   library(verification)
   library(ROCR)
   library(MASS)
   library(RColorBrewer)
   library(heatmap.plus)

   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(m[1,])
   dim(m)
   sample.names <- dataset$names

   n.models <- length(m[,1])
   temp <- strsplit(input.ds, split="/") # Extract test file name
   s <- length(temp[[1]])
   test.file.name <- temp[[1]][s]
   temp <- strsplit(test.file.name, split=".gct")
   test.file.prefix <-  temp[[1]][1]
   char.res <-  0.013 * n.models + 0.65

   # normalize scores

   if (normalize.score == T) {
     if (normalization.type == "zero.one") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
         }
     } else if (normalization.type == "z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
         }
     } else if (normalization.type == "r.z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
          }
     }         
   }

   CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   if (is.vector(cls.labels)) {
      n.phen <- 1
   } else {
      n.phen <- length(cls.labels[,1])
   }
   if (!is.na(user.colors)) {
      c.test <- user.colors
    } else {
      if (!is.null(CLS$col.phen)) {
         c.test <- CLS$col.phen
      } else {
         c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"))
      }
    }


   if (!is.null(CLS$phen.names)) {
      phen.names <- CLS$phen.names
   } else {
      phen.names <- "NA"
   }

   cls.phen.index <- unlist(cls.phen)
   cls.phen.colors <- c.test[1:length(cls.phen.index)]

   n.classes <- vector(length=n.phen, mode="numeric")
   if (n.phen == 1) {
      max.classes <- length(cls.phen)
      n.classes[1] <- max.classes
   } else {
     max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
     for (i in 1:n.phen) {
       n.classes[i] <- length(cls.phen[[i]])
     }
   }

   filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", sep="")
    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
#   windows(width=12, height=8)

   loc <- match(model, model.names)
   print(c("loc:", loc))
   s.order <- order(m[loc,], decreasing = decreasing.order)
   m2 <- m[, s.order]

   sample.names2 <- sample.names[s.order]

   if (is.vector(cls.labels)) {
      cls.labels2 <- cls.labels[s.order]
      cls.list2 <- cls.list[s.order]
   } else {
      cls.labels2 <- cls.labels[, s.order]
      cls.list2 <- cls.list[, s.order]
   }
      # Recompute cls.phen and cls.labels2 as order may have changed

     cls.phen2 <- NULL
     if (is.vector(cls.labels)) {
        classes <- unique(cls.list2)
        cls.phen2 <- classes
        cls.labels2 <- match(cls.list2, cls.phen2)
      } else {
         for (kk in 1:length(cls.list2[, 1])) {
            classes <- unique(cls.list2[kk,])
            cls.phen2[[kk]] <- classes
            cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
         }
   }


   correl <- cor(t(m2))[, loc]
   m.order <- order(correl, decreasing=T)
   correl2 <- correl[m.order]
   m2 <- m2[m.order,]
   model.names2 <- model.names[m.order]
   model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
   phen.list <- unlist(cls.phen2)
   colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
   
   MSIG.HeatMapPlot.7(V = m2, row.names = model.names2,
                      row.names2 = model.descs2, col.labels = cls.labels2, 
                      col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
                      col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = T,  
                      cmap.type = 3, char.rescale = 1,  legend=T)

   dev.off()
    
 }


OPAM.sort.projection.by.score <- function(
    input.ds,
    input.cls,
    results.dir,
    normalize.score = T,
    normalization.type = "zero.one",
    model,
    user.colors = NA,
    decreasing.order = T,
    output.dataset = NA)
  {

   library(gtools)
   library(verification)
   library(ROCR)
   library(MASS)
   library(RColorBrewer)
   library(heatmap.plus)

   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(m[1,])
   dim(m)
   sample.names <- dataset$names

   n.models <- length(m[,1])
   temp <- strsplit(input.ds, split="/") # Extract test file name
   s <- length(temp[[1]])
   test.file.name <- temp[[1]][s]
   temp <- strsplit(test.file.name, split=".gct")
   test.file.prefix <-  temp[[1]][1]
   char.res <-  0.013 * n.models + 0.65

   # normalize scores

   if (normalize.score == T) {
     if (normalization.type == "zero.one") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
         }
     } else if (normalization.type == "z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
         }
     } else if (normalization.type == "r.z.score") {
         for (i in 1:n.models) {
            m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
          }
     }         
   }

   CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   if (is.vector(cls.labels)) {
      n.phen <- 1
   } else {
      n.phen <- length(cls.labels[,1])
   }
   if (!is.na(user.colors)) {
      c.test <- user.colors
    } else {
      if (!is.null(CLS$col.phen)) {
         c.test <- CLS$col.phen
      } else {
         c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"))
      }
    }


   if (!is.null(CLS$phen.names)) {
      phen.names <- CLS$phen.names
   } else {
      phen.names <- "NA"
   }

   cls.phen.index <- unlist(cls.phen)
   cls.phen.colors <- c.test[1:length(cls.phen.index)]

   n.classes <- vector(length=n.phen, mode="numeric")
   if (n.phen == 1) {
      max.classes <- length(cls.phen)
      n.classes[1] <- max.classes
   } else {
     max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
     for (i in 1:n.phen) {
       n.classes[i] <- length(cls.phen[[i]])
     }
   }

   filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
   pdf(file=paste(filename, ".pdf", sep=""), height = 8.5, width = 11)
#   windows(width=12, height=8)

   loc <- match(model, model.names)
   print(c("loc:", loc))
   s.order <- order(m[loc,], decreasing = decreasing.order)
   m2 <- m[, s.order]

   sample.names2 <- sample.names[s.order]

   if (is.vector(cls.labels)) {
      cls.labels2 <- cls.labels[s.order]
      cls.list2 <- cls.list[s.order]
   } else {
      cls.labels2 <- cls.labels[, s.order]
      cls.list2 <- cls.list[, s.order]
   }
      # Recompute cls.phen and cls.labels2 as order may have changed

     cls.phen2 <- NULL
     if (is.vector(cls.labels)) {
        classes <- unique(cls.list2)
        cls.phen2 <- classes
        cls.labels2 <- match(cls.list2, cls.phen2)
      } else {
         for (kk in 1:length(cls.list2[, 1])) {
            classes <- unique(cls.list2[kk,])
            cls.phen2[[kk]] <- classes
            cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
         }
   }


   correl <- cor(t(m2))[, loc]
   m.order <- order(correl, decreasing=T)
   correl2 <- correl[m.order]
   m2 <- m2[m.order,]
   model.names2 <- model.names[m.order]
   model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
   phen.list <- unlist(cls.phen2)
   colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
   
   MSIG.HeatMapPlot.7(V = m2, row.names = model.names2,
                      row.names2 = model.descs2, col.labels = cls.labels2, 
                      col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
                      col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = T,  
                      cmap.type = 3, char.rescale = 1,  legend=T)

   dev.off()

   if (!is.na(output.dataset)) {
      V.GCT <- m2
      colnames(V.GCT) <- sample.names2
      row.names(V.GCT) <- model.names2
      write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename =output.dataset)  
   }
    
 }

MSIG.Subset.Dataset.2 <- function(
   input.ds,
   input.cls = NULL,
   column.subset = "ALL",    # subset of column numbers or names (or phenotypes)
   column.sel.type = "samples",  # "samples" or "phenotype"
   row.subset = "ALL",       # subset of row numbers or names
   output.ds,
   output.cls = NULL) {

# start of methodology

   print(c("Running MSIG.Subset.Dataset... on GCT file:", input.ds))
   print(c("Running MSIG.Subset.Dataset... on CLS file:", input.cls))

# Read input datasets

   dataset <- MSIG.Gct2Frame(filename = input.ds)
   m <- data.matrix(dataset$ds)
   gs.names <- dataset$row.names
   gs.descs <- dataset$descs
   sample.names <- dataset$names

# Read CLS file

   if (!is.null(input.cls)) {
      CLS <- MSIG.ReadPhenFile(file=input.cls)
      class.labels <- CLS$class.v
      class.phen <- CLS$phen
      class.list <- CLS$class.list 
   }

# Select desired column subset

   if (column.sel.type == "samples") {
      if (column.subset[1] == "ALL") {
         m2 <- m
         sample.names2 <- sample.names
         if (!is.null(input.cls)) {
            class.labels2 <- class.labels
         }
      } else {
         if (is.numeric(column.subset[1])) {
            m2 <- m[,column.subset]
            sample.names2 <- sample.names[column.subset]
            if (!is.null(input.cls)) {
              if (is.vector(class.labels)) {
                class.labels2 <- class.labels[column.subset]
              } else {
                class.labels2 <- class.labels[, column.subset]
              }
            }
         } else {
            locations <- !is.na(match(sample.names, column.subset))
            sample.names2 <- sample.names[locations]
            m2 <- m[, locations]
            if (!is.null(input.cls)) {
               if (is.vector(class.labels)) {
                  class.labels2 <- class.labels[locations]
               } else {
                  class.labels2 <- class.labels[, locations]
               }
             }
         }
      }
   } else if (column.sel.type == "phenotype") {
       locations <- !is.na(match(class.list, column.subset))
       sample.names2 <- sample.names[locations]
       m2 <- m[,locations]
       if (!is.null(input.cls)) {
          if (is.vector(class.labels)) {
             class.labels2 <- class.labels[locations]
           } else {
             class.labels2 <- class.labels[, locations]
           }
        }
   }
 
   if (row.subset[1] == "ALL") {
      m3 <- m2
      gs.names2 <- gs.names
      gs.descs2 <- gs.descs
   } else {
       locations <- !is.na(match(gs.names, row.subset))
       m3 <- m2[locations,]
       gs.names2 <- gs.names[locations]
       gs.descs2 <- gs.descs[locations]
   }

# Save datasets

   V <- data.frame(m3)
   names(V) <- sample.names2
   row.names(V) <- gs.names2
   write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)  

   if (!is.null(input.cls)) {
      write.cls.2(class.v = class.labels2, phen = class.phen, filename = output.cls) 
   }
}

OPAM.match.projection.to.pathway  <- function(
    input.ds,
    input.cls          = NA,
    results.dir,
    normalize.score    = F,
    normalization.type = "zero.one",
    pathway,
    max.n              = 10,
    user.colors        = NA,
    decreasing.order   = T,
    sort.columns       = F,
    char.rescale       = 1.25,
    cmap.type          = 3,
    row.norm           = T,
    output.dataset     = NA)
{
   library(gtools)
   library(verification)
   library(ROCR)
   library(MASS)
   library(RColorBrewer)
   library(heatmap.plus)

   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   pathway.names <- dataset$row.names
   pathway.descs <- dataset$descs
   Ns <- length(m[1,])
   dim(m)
   sample.names <- dataset$names

   n.pathways <- length(m[,1])
   temp <- strsplit(input.ds, split="/") # Extract test file name
   s <- length(temp[[1]])
   test.file.name <- temp[[1]][s]
   temp <- strsplit(test.file.name, split=".gct")
   test.file.prefix <-  temp[[1]][1]
#   char.res <-  0.013 * n.pathways + 0.65

   # normalize scores

   if (normalize.score == T) {
     if (normalization.type == "zero.one") {
         for (i in 1:n.pathways) {
            m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
         }
     } else if (normalization.type == "z.score") {
         for (i in 1:n.pathways) {
            m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
         }
     } else if (normalization.type == "r.z.score") {
         for (i in 1:n.pathways) {
            m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
          }
     }         
   }


   loc <- match(pathway, pathway.names)
   print(c("loc:", loc))
   if (sort.columns == T) {
      s.order <- order(m[loc,], decreasing = decreasing.order)
      m2 <- m[, s.order]
      sample.names2 <- sample.names[s.order]
   } else {
      m2 <- m
      sample.names2 <- sample.names
   }
   correl <- cor(t(m2))[, loc]
   m.order <- order(correl, decreasing=T)
   correl2 <- correl[m.order]
   m2 <- m2[m.order[1:max.n],]
   pathway.names2 <- pathway.names[m.order]
   pathway.descs2 <- signif(correl2, digits=3)
 
   if (input.cls == "NA") {
      cls.labels2 <- c(rep(0, 10), rep(1, length(sample.names2) - 10))
      cls.phen2 <- c(" ")
      colors.list <- c("white")
      phen.names2 <- "    "
   } else {
      CLS <- MSIG.ReadPhenFile(file = input.cls) # Read phenotype file (CLS format)
      cls.labels <- CLS$class.v
      cls.phen <- CLS$phen
      cls.list <- CLS$class.list 
      if (!is.null(CLS$phen.names)) {
         phen.names <- CLS$phen.names
      } else {
         phen.names <- "  "
      }
      if (is.vector(cls.labels)) {
         if (sort.columns == T) {
             cls.labels2 <- cls.labels[s.order]
             cls.list2 <- cls.list[s.order]
          } else {
             cls.labels2 <- cls.labels
             cls.list2 <- cls.list
          }
          n.phen <- 1
       } else {
         if (sort.columns == T) {
             cls.labels2 <- cls.labels[, s.order]
             cls.list2 <- cls.list[, s.order]
          } else {
             cls.labels2 <- cls.labels
             cls.list2 <- cls.list
          }
         n.phen <- length(cls.labels2[,1])
      }
     cls.phen2 <- NULL
     if (is.vector(cls.labels2)) {
        classes <- unique(cls.list2)
        cls.phen2 <- classes
        cls.labels2 <- match(cls.list2, cls.phen2)
      } else {
         for (kk in 1:length(cls.list2[, 1])) {
            classes <- unique(cls.list2[kk,])
            cls.phen2[[kk]] <- classes
            cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
         }
      }
      phen.names2 <- phen.names
      if (!is.na(user.colors[1])) {
         c.test <- user.colors
      } else {
         if (!is.null(CLS$col.phen)) {
            c.test <- CLS$col.phen
         } else {
            c.test <- c(brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"),
                     brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
                     brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
                     brewer.pal(n=8, name="BuGn"))
         }
       }
   }
   cls.phen.index <- unlist(cls.phen2)
   colors.list <- c.test[1:length(cls.phen.index)]

   filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ.TO.", pathway, sep="")
   pdf(file=paste(filename, ".pdf", sep=""), height = 8.5, width = 10.5)

   MSIG.HeatMapPlot.7(V = m2, row.names = pathway.names2,
                      row.names2 = pathway.descs2, col.labels = cls.labels2, 
                      col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
                      col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = row.norm,  
                      cmap.type = cmap.type, char.rescale = char.rescale,  legend=T)
   dev.off()

   if (!is.na(output.dataset)) {
      V.GCT <- m2
      colnames(V.GCT) <- sample.names2
      row.names(V.GCT) <- pathway.names2
      write.gct(gct.data.frame = V.GCT, descs = pathway.descs2, filename =output.dataset)  
   }
    
}
