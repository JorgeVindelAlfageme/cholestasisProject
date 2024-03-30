## 1) checkPackage:

# It attaches or installs, if it is not already on the computer, an R library.

checkPackage <- function(packageLabel, # String of a package name to be
                         # attached or installed.
                         repository = "CRAN" # String of the R repository the
                         # previous package can be found in and installed.
                         ){
  
  if(!(packageLabel %in% rownames(installed.packages()))){
    # checks if the package is among the computer installed packages. In case
    # it is not installed, this code chunk is executed.
    
    if(repository == "CRAN"){ # If the repository where the package is in is
      # CRAN, this piece of code is run:
      catMe <- paste("Installing ", packageLabel, " from CRAN...\n", sep = "")
      cat(catMe)
      install.packages(packageLabel) # CRAN package installation is run.
    }
    if(repository == "Bioconductor"){ # If the repository where the package is
      # in is Bioconductor, this piece of code is run:
      catMe <- paste("Installing ", packageLabel, " from Bioconductor...\n", sep = "")
      cat(catMe)
      BiocManager::install(packageLabel) # Bioconductor package installation is
      # run.
    }
    
  }else{ # If the package name can be found among the packages that are already
    # installed on the computer, a message indicating so is emitted.
    
    catMe <- paste(packageLabel, " is already installed.\n", sep = "")
    cat(catMe)
    
  }
  
  packageCode <- paste("package::", packageLabel) # Installed package tag is
  # recorded.
  
  if(!(packageCode %in% search())){ # Package tag is scanned among already
    # attached packages. If the package has not been attached yet, this first
    # piece of code is executed.
    
    library(packageLabel, character.only = TRUE) # attaches indicated package.
    catMe <- paste(packageLabel, " successfully attached.\n", sep = "")
    cat(catMe)
    
  }else{ # If the package is already attached, a message indicating is already
    # loaded is printed.
    
    catMe <- paste(packageLabel, " was already attached.\n", sep = "")
    cat(catMe)
    
  }
  
}

## 2) geneIDextraction:

# This function was built to get genes IDs from a protein description generated
# from Proteome Discoverer (TM) software, which at the same time was obtained
# from UniProt database. It makes use of a regular expression to extract the
# gene name of a given protein.

geneIDextraction <- function(geneDescription # This function requires a protein
                             # description as starting point in which a
                             # particular character code can be found. From
                             # this part of the description, the gene name of
                             # the protein can easily be extracted. Input can
                             # be a character vector to process multiple
                             # descriptions.
                             ){
  
  genePattern <- "GN=(.*?) " # encrypts wanted gene name pattern as a regular
  # expression.
  geneID <- regmatches(geneDescription, regexpr(genePattern, geneDescription))
  # finds regular expression pattern from a character vector and from string
  # indexing to retrieve all substring matches and saves them as another
  # character vector.
  geneID <- gsub(pattern = "GN=", replacement = "", x = geneID) # erases part
  # of the regular expression used to locate wanted substrings which is not
  # necessary anymore.
  geneID <- gsub(pattern = " ", replacement = "", x = geneID) # deletes another
  # part of the substring, an unwanted blank space.
  
  # Nevertheless, some protein descriptions lack their corresponding gene name,
  # which produces a shortening in retrieved gene names character vector. To
  # make sure the produced gene name character vector has the same length as
  # the original given protein description compilation, a loop iterating over
  # the entire protein description compilation will be run. Thanks to this
  # loop, absent names will be filled with a string corresponding an empty name.
  # Proteins lacking their gene names in the description will be given an empty
  # string as name, but at least it will not be completely nonexistent.
  
  realGeneID <- c() # empty vector that will hold gene names and which length
  # will equal protein description character vector given as an input.
  j <- 1 # creates an index which will be updated after each loop iteration,
  # starting at one.
  for(i in 1:length(geneDescription)){ # This loop iterates over the original
    # protein descriptions vector.
    istr <- geneDescription[i] # retrieves one exclusive whole protein
    # description.
    if(grepl(geneID[j], istr, fixed = TRUE)){ # If the extracted gene name
      # following the gene vector order can be found in the description, this
      # part of the code is executed.
      realGeneID <- append(x = realGeneID, values = geneID[j]) # Current gene
      # name is added to the previous empty vector.
      j <- j + 1 # Gene name vector index is increased by one unit.
    }else{ # If the currently selected gene name cannot be found in the protein
      # description, it will be supposed that the gene name is absent in that
      # description.
      realGeneID <- append(x = realGeneID, values = "") # A blank string will
      # be added to the output vector just so the original protein descriptions
      # vector and the final gene names vector have the same length.
    }
  }
  
  return(realGeneID) # outputs genes names vector.
}

## 3) PCACal:

# A simple function that automates Principal Components Analysis calculation
# based on "prcomp" built-in function.

PCACal <- function(scaledData, # Data matrix with its samples arranged by rows
                   # and its features disposed by columns.
                   wanna.scale = FALSE # "prcomp" implicitly scales data before
                   # executing the PCA by default. This argument is passed
                   # directly to "prcomp" scaling argument, to control whether
                   # or not data should be normalized.
                   ){
  
  PCA <- prcomp(x = t(na.omit(scaledData)), scale. = wanna.scale) # executes
  # PCA.
  return(PCA) # returns PCA.
  
}

## 4) elbowM:

# It executes the elbow method based on a series of values contained in a
# vector. The elbow method is a heuristic method usually used to determine
# the optimal number of clusters based on magnitudes such as sum of squared
# distances when attending spatial features in clustering procedures. It is
# also used on PCA variance retention to decide the best PCs to be plotted. 

elbowM <- function(vec, # Vector containing a series of values. This set of
                   # values should be characterized by that the first one is
                   # the greatest of them, and values are progressively
                   # smaller.
                   namedVec = FALSE # A boolean argument to indicate if
                   # previous vector is a named vector.
                   )
  {
  
  ratesVec <- c() # Vector that accumulates rates by dividing two consecutive
  # vector differences. Differences are calculated by substracting two
  # consecutive numbers, particularly the smallest one to the greatest one.
  # Rates are calculated by dividing differences that have been calculated in
  # order, starting from the greatest values to the smallest ones.
  diffsVec <- c() # Vector that will gather vector differences.
  len <- length(vec) # Length of original values vector.
  for(i in 2:len){ # A loop that iterates over indices, starting at 2 and
    # finishing at a number that equals the length of the original values
    # vector.
    if(i <= (len - 1)){ # If the current loop index is less than the length of
      # the original values vector minus one, this part of the code is run.
      rate1 <- unname(vec[i-1]) - unname(vec[i]) # calculates difference
      # between previous vector element and currently selected vector element.
      # The variable name used is because this result will be used to calculate
      # a rate between two differences.
      rate2 <- unname(vec[i]) - unname(vec[i+1]) # calculates difference
      # between currently selected vector element and subsequent vector element.
      rate <- round(rate1/rate2, 5) # calculates rate between last two obtained
      # differences. Also, result is rounded to 5 decimal places.
      ratesVec <- append(x = ratesVec, values = rate) # computed rate is added
      # to rates vector.
    }
    
    diff <- vec[i-1] - vec[i] # computes difference between two consecutive
    # vector values. It is equivalent to the result of "rate1" variable. This
    # variable differs in its name and does not constitute a named vector.
    diffsVec <- append(x = diffsVec, values = diff) # appends difference
    # result to differences vector.
    
  }
  
  ratesVec <- c(NA, ratesVec, NA) # adds two "NA" ("not available") values to
  # rates vector: one at its beginning and one at its end. This way, the length
  # of this vector and the original values vector is the same. It was thought
  # to represent all results from this function as a data frame at the end, so
  # making vectors having the same length makes merging them into a single data
  # frame a logical approach.
  diffsVec <- c(NA, diffsVec) # lengthens differences vector by adding one "NA"
  # value at its beginning.
  dataF <- rbind(vec, ratesVec, diffsVec) # creates data frame from original
  # values vector, differences vector and rates vector.
  rownames(dataF) <- c("Original vector",
                       "Continuous difference rate",
                       "Continuous difference")
  # names previous data frame rows.
  
  if(namedVec == TRUE){ # If it was indicated that the original values vector
    # was named, this part of the conditional is run.
    colnames(dataF) <- names(vec) # names data frames columns after original
    # values vector names.
  }else{ # If it was indicated that the original values vector was not named,
    # this code chunk is run:
    colnames(dataF) <- 1:len # just names columns data frame using indices from
    # one to the total number of columns it has.
  }
  print(dataF) # prints data frame on command prompt.
  orderedVec <- diffsVec[-1][order(diffsVec[-1], decreasing = TRUE)]
  # creates numerically ordered version of differences vector.
  if(all(as.numeric(diffsVec[-1]) == as.numeric(orderedVec))) # If the
    # numerically ordered differences vector fully coincide with its original
    # form, it is recommended to plot those components which are related to
    # maximization of differences rates values.
    {
    return(paste0("Recommended plotting components: ",
                  as.character(colnames(dataF)[which.max(ratesVec)-1]),
                  ", ",
                  as.character(colnames(dataF)[which.max(ratesVec)]),
                  sep = ""))
  }else{ # If the ordering between vector values does not coincide,
    # maximization of differences values will be considered to determine
    # heuristically optimal number of clusters, selected components, etc.
    return(paste0("Recommended plotting components: ",
                  as.character(colnames(dataF)[which.max(diffsVec)-1]),
                  ", ",
                  as.character(colnames(dataF)[which.max(diffsVec)]),
                  sep = ""))
  }
}

## 5) PCAPlot:

# It plots PCA results from previous "PCACal" function into a 2D scatter graph,
# with each point representing a sample from data.

PCAPlot <- function(PCAObject, # Variable returned by "prcomp" function.
                    is.chosen = FALSE, # refers to the Principal Components (PC)
                    # that are wanted to be plotted. If "TRUE", PCs have to be
                    # chosen by their indices in the next arguments. If "FALSE",
                    # this function will decide which PCs to plot based on
                    # the variance retained as a percentage or proportion by
                    # each PC and the elbow method results upon variance values.
                    firstPC, # If "is.chosen = TRUE", this argument serves to
                    # indicate chosen PC results to be represented across the
                    # horizontal axis.
                    secPC, # If "is.chosen = TRUE", this argument serves to
                    # indicate chosen PC results to be represented across the
                    # vertical axis.
                    groups, # Character vector with the name of the classes
                    # that are also contained as substrings in samples names.
                    # Substring belonging analysis will be used to attribute
                    # a color to each sample from a class, and all this is 
                    # based on samples substring identification.
                    colVec, # A character vector with the same length as the
                    # previous classes vector to associate each class to a
                    # color. Therefore, both classes vector and this colors
                    # vector have to be ordered so each class matches its color.
                    legPos = "topright", # This plot will produce a legend.
                    # This argument controls legend position in the plot frame.
                    main2 = "default", # This argument is used to create a plot
                    # title. If it equals "default", a default title is
                    # generated, indicating the scatter plot is a PCA and which
                    # PCs are being plotted. Any other title can be printed
                    # just by writing the desired title as a string.
                    cexLeg = 1, # Factor that controls legend size.
                    cexPoints = 1 # Factor that controls plot points size.
                    ){
  if(is.chosen == FALSE){ # If it is wanted that the function considers what
    # PCs are wished to be plotted based on the elbow method results on
    # variance retained per principal component, it should be indicated
    # through the "is.chosen" argument as a boolean value and this part of
    # the code will be executed.
    
    elbow <- elbowM(summary(PCAObject)$importance[2,], namedVec = FALSE)
    # computes elbow method results based on previously shown user-defined
    # function.
    
    firstPC <- as.numeric(substr(x = elbow,
                                 start = nchar(elbow)-3,
                                 stop = nchar(elbow)-3)) # retrieves the first
    # recommended principal component to be plotted by the returned string as
    # output the elbow method defined function generates.
    secPC <- as.numeric(substr(x = elbow,
                               start = nchar(elbow),
                               stop = nchar(elbow))) # retrieves the second
    # recommended principal component by the elbow method function. 
  }
  
  if(main2 == "default"){ # If no specific plot title was written down, a
    # default one is produced.
    
    main2 <- paste0("PC",
                    as.character(firstPC),
                    " vs. PC",
                    as.character(secPC),
                    sep = "") # Default plot title will be based on chosen
    # principal components, which will be identified in it using their indices.
    
  }
  
  plot(x = PCAObject$x[,firstPC], # plots first chosen Principal Component
       # values across horizontal axis.
       y = PCAObject$x[,secPC], # plots second chosen Principal Component
       # values across vertical axis.
       xlim = range(pretty(c(PCAObject$x[,firstPC]))), # calculates maximum and
       # minimum value from a set of values that allow aesthetic graphical
       # representation generated from plot horizontal axis values.
       ylim = range(pretty(c(PCAObject$x[,secPC]))), # calculates maximum and
       # minimum value of a set of values that allow aesthetic graphical
       # representation generated from plot vertical axis values.
       main = main2, # Plot title.
       xlab = paste0("PC", as.character(firstPC), " (",
                     as.character(round(100*summary(PCAObject)$importance[2,firstPC],1)),
                     "%)",
                     sep = ""), # Horizontal axis name label. What is being
       # added as part of the label is the variance retained by the plotted
       # Principal Component represented as a percentage of the full data set.
       ylab = paste0("PC", as.character(secPC), " (",
                     as.character(round(100*summary(PCAObject)$importance[2,secPC],1)),
                     "%)",
                     sep = ""), # Vertical axis name label. What is being
       # added as part of the label is the variance retained by the plotted
       # Principal Component represented as a percentage of the full data set.
       cex = cexPoints # Factor that controls plot points size.
       ) # plots PCA results as a scatter plot.
  
  for(i in 1:length(groups)){ # A loop that iterates over the unique class tags
    # character vector. It selects one class tag at a time to later search for
    # the data set samples that can be grouped into that class. By doing so,
    # samples are colored following a distinctive pattern.
    
    idx <- which(grepl(pattern = groups[i], x = rownames(PCAObject$x), fixed = TRUE))
    # Samples from the different classes are identified as indices from the
    # column names of the data set. Each subset of samples, depending of its
    # class, will be specifically colored.
    
    points(x = PCAObject$x[idx, firstPC], # Horizontal axis values for selected
           # class samples.
           y = PCAObject$x[idx, secPC], # Vertical axis values for selected
           # class samples.
           pch = 21, # Type of point that is wanted to be used.
           col = "black", # Color of points perimeter.
           bg = colVec[i], # Color chosen to represent this subset of
           # samples points considering their class.
           cex = cexPoints # Factor that controls points size.
           )
    
  }
  
  legend(x = legPos, # places the legend in a predetermined quadrant position.
         legend = groups, # Names of the classes that will be used as the
         # legend text.
         pch = 21, # Type of marker that is going to be used as samples points. 
         col = "black", # Color for makers margin.
         pt.bg = colVec, # Original colors vector matching their corresponding
         # classes tags.
         cex = cexLeg # Factor that controls legend size.
         )
  
}

## 6) fortuna:

# "Fortuna" represents a module used for feature selection, which would be
# specially useful when the number of features is gargantuan. It is based on the
# feature importance measurement random forest provides. Iteratively, the least
# important features will be deleted from the analysed data set. Furthermore,
# a comparison between the kept features in the algorithm and their randomized
# version will be carried out at the same time. Also those features that are
# kept in the data set for successive rounds are considered to be more important
# than the ones that have been already eliminated. Once there is only one sole
# feature left in the data set, the entire set will be regenerated and the
# process of feature elimination and comparison will be repeated as many times
# as the user decides. From this information, a final feature importance value
# and a p-value out of a binomial test will be output for each feature,
# establishing a feature ranking.

# Fortuna was envisioned after Diaz-Uriarte experiments (2006) and I did not
# know the existence of Boruta (Kursa, 2010), which is also based on random
# forest and recursive feature elimination. Boruta works incredibly faster
# than Fortuna, and statistical results rely on multiple chained tests.
# Nevertheless, Fortuna brings data set regeneration and the repetition of the
# process could be a included step to increase results reliability.

# Fortuna also brings innovative approaches for selecting features by altering
# the way it chooses them by changing the number of features eliminated per
# round. It allows the user to modify the number of features eliminated per
# round by selecting one out of two different modes. One of them is based on
# indicating the proportion of features that will be retained after each
# elimination round. The second one is based on indicating the number of
# elimination rounds that will be in total starting from the entire data set
# until there is only one feature left. During every iteration, importance
# measurement and feature randomization comparison and feature retention
# assessment will be carried out.

# Furthermore, Fortuna brings an innovative way of grouping samples and
# obtaining importance results, based on how samples are put into groups.
# One way of analyzing importance is just by considering each different class
# at the same time. Another way of computing the importance analysis is by
# opposing each class to the rest of samples under a unifying class tag.
# That way, the number of importance analyses and results will equal the
# number of classes in the study. Emitting multiple binary classification
# analyses can stress features that clearly can distinguish classes under
# certain conditions that could not be found by other approaches.

# Finally, after obtaining the importance results, this algorithm will
# elaborate an eventual random forest classifier with the features that have
# been considered important to classify the samples. This last classifier (or
# set of classifiers if the classification approach has been selected to be
# binary) will have random forest proximity measurements and also feature
# importance measurements to extract further knowledge from it.

fortuna <- function(dataset, # Matrix or data frame having their rows
                    # representing samples, and their columns representing
                    # features.
                    groups, # Vector being a factor object indicating the
                    # classes the samples from the previous matrix or data
                    # frame belong to.
                    treeNumber = 1000, # Number of trees each random forest
                    # will have.
                    retentionByProp = FALSE, # indicates what kind of feature
                    # elimination is going to be produced. If "TRUE", feature
                    # elimination would consider in each elimination round
                    # a proportion of the most important features that will
                    # be kept to be retained in the algorithm. If "FALSE",
                    # the user will decide for how many round features will be
                    # erased until there is only one left by making use of
                    # another argument. Either if this argument equals TRUE or
                    # FALSE, there will be another argument to select the
                    # features retention proportion or the number of elimination
                    # rounds.
                    retentionProportion = 0.8, # Proportion of features
                    # retained in the algorithm after each elimination round.
                    retentionLoops = 30, # Number of elimination loops that
                    # will be run until there is only one feature left starting
                    # from the entire feature set.
                    chosenIterationsNumber = 30, # Number of times the entire
                    # data set will be regenerated.
                    isBinaryApproach = FALSE, # This argument tries to provide
                    # with another possible feature selection approach by
                    # reordering samples between different and compatible
                    # groups. If "isBinaryApproach = FALSE", the entire
                    # algorithm process will be only executed one, comparing
                    # all classes at the same time. The most important features
                    # will be those that are able to distinguish between all
                    # classes simultaneously. By contrast, if it equals "FALSE",
                    # each group will be opposed to the rest of the samples
                    # grouped under the same unifying class tags. What this will
                    # will try to do is to find those features that characterize
                    # each class with respect to the rest. At the end, there
                    # will be as many sets of results as classes are in the
                    # data set.
                    chosenFeatsClassMethod = "elbow", # This argument is used to
                    # indicate which criterion has to be used to compile a set
                    # of features to be used in the generation of a final
                    # random forest classifier. The selected features will
                    # always be among the most important. This argument can
                    # take 3 different character values: "k", "elbow", "p".
                    # "k" just stands for a set of the "k" most important
                    # features. "elbow" selects those most important features
                    # based on the elbow method results after ordering all
                    # importance values, from greatest to lowest. "p" stands
                    # for those significantly important features to be selected.
                    chosenPValue = 0.05, # Significance threshold for feature
                    # importance selection based on p-values.
                    pAdjMethod = "fdr", # P-value adjustment method. Possible
                    # argument values coincide with those of "p.adjust"
                    # built-in function.
                    kFeats = 20, # Number of most important features to be
                    # chosen for classifier building if "chosenFeatsClassMethod"
                    # equals "k".
                    isElbowPlot = TRUE # can take "TRUE" of "FALSE" values.
                    # Feature importance values are decreasingly ordered and
                    # the elbow method is applied on them to distinguish two
                    # parts in the linear plot. Using a color code, important
                    # features after elbow method results are highlighted.
                    # This plot helps visualizing the proportion of features
                    # from a data set that seem relatively important over the
                    # rest.
                    ){
  
  if(!("package:randomForest" %in% base::search())){
    library(randomForest)
  } # loads "randomForest" library if it has not been loaded yet.
  
  if(!(retentionByProp)){ # If feature elimination by data set feature
    # proportion retention has not been selected as the feature elimination
    # method, this part of the code will be run.
    
    # Basically, from a given proportion and considering the number of total
    # features, the number of elimination rounds can also be calculated. This
    # number can just be introduced in a loop as the number of iterations.
    # Eventually, any feature elimination approach conceded by this function
    # is related to computing the feature elimination proportion.
    
    howManyLoops <- function(nfeat, wantedLoops){ # calculates the proportion
      # that is necessary to reduce a number to one by iteratively
      # substracting it starting from a total initial number.
      prob <- seq(1 - 0.001, 0.001, by = -0.001) # A vector of numbers between
      # almost one and zero to approximate the number of rounds to finally
      # iteratively reduce the original feature number to one.
      iterationsVector <- c() # Empty vector that will hold the number of
      # elimination iterations each proportion retention results in. After
      # executing this code, a set of proportions will be adequate to be used
      # as feature retention proportions. The obtained numbers in this vector
      # have to match with the one chosen by the user at the function
      # execution.
      for(i in 1:length(prob)){ # Loop that iterates over all the elements of
        # generated proportions.
        remaining <- c() # Empty vector that will hold the reduced number that
        # is left after multiplying it by the currently selected proportion
        # after each inner loop iteration.
        # After the original number has been converted into one, the length of
        # this vector is calculated and after using the correct proportion its
        # length will equal the number of desired elimination rounds.
        # It empties at each outer loop iteration.
        remember <- nfeat # This value is updated after each inner loop
        # iteration. It is consecutively reduced by the product of it and the
        # currently selected proportion value. It updates its value to the
        # number of total features at the beginning of each outer loop
        # iteration. 
        while(remember > 1){ # This loop will be executed while the
          # consecutively reduced number by a product of it and the selected
          # proportion is greater than one.
          remember <- floor(remember * prob[i]) # It calculates the product of
          # the remaining iteratively reduced number and the proportion that
          # consecutively makes it smaller. Plus, the result is always rounded
          # to the smallest integer.
          remaining <- append(x = remaining, values = remember) # incorporates
          # recently calculated product to a vector.
        }
        iterationsVector <- append(x = iterationsVector, values = length(remaining))
        # After the inner loop has finished, the length of the vector containing
        # all the results of the proportionately diminished total number of
        # features is computed and added to another vector.
      }
      
      return(prob[iterationsVector == wantedLoops]) # Only the proportions
      # which implied the wanted number of reducing rounds will be output.
      # Any number from this output can be used to convert the number of total
      # features into one.
      
    }
    
    retentionProportion <- howManyLoops(nfeat = ncol(dataset),
                                        wantedLoops = retentionLoops)[1]
    # calculates the retention proportion a number is linked to to eventually
    # equal one if it is multiplied by the proportion a wanted number of times.
    
    # It is possible that the theoretical wanted number of elimination rounds
    # exceeds the number of features. In that case, previous function will
    # output a "NA" value. In that case, the elimination proportion will equal
    # 0.99 and the number of elimination rounds will be equal to the number of
    # features in the given data set by default. What this will produce is that
    # in each elimination round only the least important feature will be
    # eliminated.
    
    if(is.na(retentionProportion)){ # If the retention proportion to reduce a
      # number to one by multiplying it for a desired number of rounds does not
      # theoretically exist, a NA value is generated and this part of the code
      # is executed.
      cat(paste("The number of features in the chosen data set is relatively ",
                "small. A single feature will be eliminated by step instead.\n",
                sep = ""))
      retentionProportion <- 0.99 # The retention proportion for the iterative
      # reduction process is scripted to be 0.99.
      chosenIterationsNumber <- ncol(dataset) # The number of elimination
      # rounds will be equal to the number of features in the currently
      # analyzed data set.
    }
    
  }
  
  # The next function is part of the module to obtain results to finally apply
  # a binomial test to get a value which will statistically represent if a
  # feature is important compared to the rest of features in the data set. It
  # will be used to generate input values to a binomial test by generating
  # randomized versions of features and comparing them to their original form.
  # Only features which importance value is greater than their chaotic version
  # will gain statistical significance power after each iterative step in the
  # analysis.
  
  fortunaRandomizer <- function(dat, # Input data set.
                                wannaCBind = TRUE # Boolean value to indicate
                                # if the result should be a merged array
                                # considering both the original one and its
                                # randomized version or not.
                                ){
    
    rndmArr <- apply(X = dat, # Input data.
                     MARGIN = 2, # The given function is used on columns.
                     FUN = function(x) sample(x) # randomizes vector elements
                     # positions.
                     )
    # generates a randomized version of original data by shuffling values
    # across rows (i. e., values are jumbled by columns). To do so, "apply"
    # function is used and its "MARGIN" argument equals 2.
    
    colnames(rndmArr) <- sprintf("rndm.%s", colnames(dat)) # names randomized
    # array columns after original column names and by adding the "rndm."
    # substring in front every column name.
    
    # "wannaCBind" is an argument that can take either "TRUE" or "FALSE" values.
    # If it equals TRUE, both the original data matrix and its generated
    # randomized version will be collapsed into a single data array. This is
    # intended because the entire merged data set can now be given to
    # random forest algorithm as an input. From this set, the results of
    # importance of each feature and its randomized version can be compared.
    # From this comparison, a value that will be used in a binomial test can be
    # obtained. If the randomized version has an importance lesser than the
    # importance of its original feature, a 1 (which will also be referred as
    # a "hit") is conceded to that feature. In contrast, if the randomized
    # version of a feature is greater than its original counterpart, a 0 is
    # conceded to that feature. Those 0s ans 1s are accumulated over multiple
    # rounds and, after a substantial number of values have been generated, a
    # binomial test can be executed to distinguish between those features that
    # have more probabilities of being important from the rest.
    
    if(wannaCBind){ # If it is wanted that the randomized array and the
      # original array are merged, based on the value given to one of the
      # function arguments, this conditional will be executed.
      
      rndmArr <- cbind(dat, rndmArr) # binds the original data array columns
      # and the randomized data array columns.
      
    }
    
    return(rndmArr) # returns final data array.
    
  }
  
  # The next function iteratively eliminates the least important features by
  # rounds, accumulates binomial test input data and results from the
  # randomized features version comparison and the most important features
  # prevalence in the algorithm, regenerates the entire data set and runs the
  # iterative feature elimination process as many times as indicated by the
  # user and outputs a matrix with all results in the form of a matrix.
  
  rrfIterator <- function(innerDataSet, # Input training data set.
                          innerGroups, # Input training samples classes vector.
                          innerTreeNumber, # Number of trees each random forest
                          # will have in the algorithm.
                          innerRetentionProportion, # Proportion of retained
                          # features to eliminate those least important ones
                          # per algorithm elimination round.
                          innerRetentionLoops, # Number of loops during which
                          # features will be eliminated to preserve those that
                          # are more important.
                          innerIsElbowPlot # Boolean value to determine if a
                          # plot representing ordered feature importance values
                          # is going to be plotted ("TRUE") or not ("FALSE").
                          )
  {
    
    importanceMatrix <- matrix(data = 0, # Unique data value to fill the entire
                               # matrix shell.
                               nrow = ncol(innerDataSet), # Number of rows the
                               # matrix shell will have. It equals the number
                               # of features that the training data set has.
                               ncol = 8 # Predetermined number of columns the
                               # matrix shell will have. Each column will be
                               # destined to store a particular portion of
                               # information based on the importance of
                               # features.
                               )
    rownames(importanceMatrix) <- colnames(innerDataSet) # The row of results
    # matrix will be named after data set features.
    colnames(importanceMatrix) <- c("feature.index", # This column inserts an
                                    # index to be able to name a feature using
                                    # an alternative code as an ID.
                                    "importance", # Raw importance value
                                    # obtained from adding all importance
                                    # values from all random forest executions
                                    # while a particular feature was still
                                    # contemplated between those retained in
                                    # the algorithm.
                                    "normalized.importance", # Normalized
                                    # importance values obtained from using raw
                                    # importance values in a rate. All positive
                                    # values from the algorithm and all negative
                                    # ones were separately considered in two
                                    # different rate calculations.
                                    "standard.cumulative", # Cumulative sum of
                                    # normalized importance values starting
                                    # with the most important features and
                                    # ending with the least important ones.
                                    "binomial.pvalue", # P-value results from
                                    # binomial test.
                                    "adj.pvalue", # Adjusted p-values.
                                    "y0.elbow.method", # A recorded value
                                    # based on a linear regression to compute
                                    # which features accumulate more importance
                                    # based on the heuristic elbow method.
                                    "is.elbow.method" # A boolean vector which
                                    # shows what features are considered to be
                                    # heuristically more important based on
                                    # the elbow method results applied onto
                                    # feature importance measurements.
                                    ) # Results matrix columns names. Each
    # column will be destined to store a different part of the results
    # information. 
    importanceMatrix[, "feature.index"] <- 1:nrow(importanceMatrix)
    # An ID based on features original order in the training data set is
    # attributed to each feature.
    
    innerChosenIterationsNumber <- 0 # Number of looped iterations that will be
    # needed to make the number of original features to be reduced to one by
    # eliminating a given proportion of the least important features in each
    # round. By considering the feature retention proportion, the number of
    # loops during which features will be eliminated can be computed. This
    # number can be used to generate a matrix that will be destined to store
    # the information retrieved from the randomization and the algorithm
    # feature preservation analysis to eventually calculate the binomial test
    # results.
    remember <- ncol(innerDataSet) # This number will be updated after each
    # iteration in the next loop to compute the number of iterations during
    # which the total number of features will be reduced to one feature. This
    # number of rounds is computed by considering the proportion of features
    # that will be retained in each feature elimination round.
    while(remember > 1){ # As long as the number of represented retained
      # features is greater than one, this loop will continue to be executed. 
      remember <- floor(remember * innerRetentionProportion) # Number of
      # features that will be retained in the algorithm after eliminating those
      # less important.
      innerChosenIterationsNumber <- innerChosenIterationsNumber + 1
      # The iteration number is increased in one loop.
    }
    
    hitArr <- matrix(data = NA, # "NA" constant value to fill the entire
                     # results matrix at first.
                     nrow = innerChosenIterationsNumber * 2 * innerRetentionLoops,
                     # The number of rows has to equal the hypothetical maximum
                     # value of results a feature can store due to its maximum
                     # theoretical importance. That number has to equal the
                     # number of elimination rounds from having the entire data
                     # set until there is only one feature left times the number
                     # of times the entire data set will be regenerated times 2.
                     # The first product has to be doubled because for every
                     # feature that prevails in the algorithm due to its
                     # importance in each round, there will be 2 values that
                     # will be associated to it. One has to do with the
                     # comparison between the original feature importance and
                     # its randomized version importance. The other has to do
                     # with if the analysed feature has prevailed during the
                     # elimination rounds until the current one or if it has
                     # been already eliminated.
                     ncol = ncol(dataset), # The number of columns of this
                     # array equals the number of total features in the
                     # analysis.
                     dimnames = list(NULL, colnames(dataset)) # names columns
                     # after feature names.
                     ) # This array will contain those results that will allow
    # the binomial test p-value calculation. The number of columns will equal
    # the number of original features. Each column will contain information
    # related to one feature at a time. Next, the number of rows will equal
    # the number of elimination loops for the entire data set times the number
    # of times the entire data set will be regenerated times 2. We have to
    # consider that there will be as many rounds of the data set reduction to
    # reach just one feature as it is desired by the user. The entire process
    # of looped feature elimination is repeated because the repetition of this
    # process guarantees more robust results. The maximum possible value for
    # one hypothetical most important feature would imply that for all
    # elimination rounds across all elimination processes after each data set
    # regeneration it has been preserved in the algorithm as the most important
    # feature with no other feature to compare it to. The maximum number of
    # iterations has to be doubled because, at each iteration, there will be
    # two values associated to the binomial test results that will be stored
    # for each feature if the feature has prevailed in the algorithm. These
    # two values have to do with: 1) the feature generation by randomization
    # and 2) feature prevailing in the algorithm.
    
    iterations <- 1 # A number constituting the times the entire data set has
    # been regenerated. This number is updated by increasing in one unit after
    # each outer loop.
    while(iterations <= innerRetentionLoops){ # As long as the number of
      # input data set regeneration rounds is not exceeded, this loop will be
      # still executed.
      cat(paste("Iteration for recursive feature elimination: ",
                as.character(iterations),
                ".\n",
                sep = "")
          ) # prints the current outer loop iteration it is being run.
      tempDataSet <- innerDataSet # creates a copy of the entire training data
      # set. 
      
      retentionLoop <- 1 # represents the number of elimination rounds that
      # have occurred. It will be increased by one unit after each elimination
      # round, and it will stop increasing when there is only one feature left
      # in the data set.
      
      allFeatsEliminated <- c() # Character vector that accumulates the
      # eliminated features since the successive feature elimination rounds
      # start. It resets its content when the entire data set is regenerated.
      
      while(!(is.null(ncol(tempDataSet)))){ # This loop will be executed until
        # there is only one feature left, which implies that the training data
        # set loses its matrix or data frame class object nature and becomes a
        # numeric object.
        
        rndm_tempDataSet <- fortunaRandomizer(dat = tempDataSet, # starting
                                              # data set.
                                              wannaCBind = TRUE # indicates
                                              # that the output matrix has to
                                              # be a merged array using the
                                              # original features and their
                                              # randomized versions by columns.
                                              )
        # returns the original data set merged with randomized versions of all
        # features from the original one created by the randomization function
        # described earlier.
        
        rndm_tempRF <- randomForest(x = rndm_tempDataSet, # Merged data set
                                    # with the original and randomized features.
                                    y = innerGroups, # samples classes as a
                                    # factor class object.
                                    importance = TRUE, # It indicates to output
                                    # the feature importance measurement as
                                    # part of the result.
                                    ntree = innerTreeNumber # Number of trees
                                    # in random forest.
                                    )
        # runs random forest using the merged data set by combining the
        # original features and their randomized version.
        
        rndm_Imp <- rndm_tempRF[["importance"]][,"MeanDecreaseAccuracy"]
        # retrieves importance measurement from the random forest object of
        # all features, both the original ones and the randomized ones.
        
        tempFeats <- colnames(tempDataSet) # retrieves the name of those
        # original features that prevail in the importance measurement
        # algorithm.
        
        boolRndm <- grepl(pattern = "rndm", # Substring that has to be found in
                          # features names and indicates if the feature is a
                          # randomized feature.
                          x = names(rndm_Imp), # Vector of features names
                          # obtained from the features importance vector.
                          fixed = TRUE # indicates a substring has to be
                          # searched as it is, without considering regular
                          # expression characters.
                          ) # boolean vector indicating which features in the
        # merged data set are the randomized ones based on their given name.
        
        for(i in 1:length(tempFeats)){ # A loop that executes as many iterations
          # as features are left in the data set.
          
          ifeat <- tempFeats[i] # orderly retrieves the original name of a
          # feature.
          iboolBoth <- grepl(pattern = ifeat, # Feature name.
                             x = names(rndm_Imp), # Features names coming from
                             # importance measurement output.
                             fixed = TRUE # indicates a substring has to be
                             # searched as it is, without considering regular
                             # expression characters.
                             ) # Boolean vector indicating the position of
          # a given feature and its randomized version.
          iboolOriginal <- !(boolRndm) & iboolBoth # Intersection between
          # two boolean vectors that indicates the position of the original
          # feature in the importance results.
          iboolRndm <- boolRndm & iboolBoth # Intersection between
          # two boolean vectors that indicates the position of the randomized
          # feature version in the importance results.
          ihitOrNot <- as.numeric(rndm_Imp[iboolOriginal] > rndm_Imp[iboolRndm])
          # compares the importance value of the selected feature and its
          # randomized version and converts a boolean result into a numeric
          # result. Greater importance attributed to the original feature
          # results in a 1 (which is call a "hit"), and the opposite
          # results in a 0.
          hitArr[which(is.na(hitArr[,ifeat]))[1],ifeat] <- ihitOrNot
          # includes last calculated result (0 or 1) in the hits shell array
          # created before, that accumulates these results for all features.
        }
        
        tempRF <- randomForest(x = tempDataSet, # Data set having only original
                               # versions of features that have prevailed in
                               # the algorithm due to their relatively higher
                               # importance.
                               y = innerGroups, # Vector indicating samples
                               # classes.
                               importance = TRUE, # indicates that it is
                               # wanted to receive an importance measurement
                               # as output.
                               ntree = innerTreeNumber # Number of trees in
                               # random forest.
                               ) # runs a second random forest, only
        # considering those features that are original.
        
        tempImportance <- tempRF$importance # Importance output is isolated
        # from the rest of the last random forest output.
        
        vecImp <- tempImportance[, "MeanDecreaseAccuracy"] # retrieves features
        # importance values with no randomized features involved.
        
        importanceMatrix[rownames(tempImportance),"importance"] <-
          importanceMatrix[rownames(tempImportance),"importance"] + vecImp
        # It incorporates the importance results from this iteration to the
        # importance results matrix shell created before, as a summation
        # of importance over loops.
        
        newOrder <- order(vecImp, decreasing = TRUE) # computes features order
        # considering their importance value which has just been calculated.  
        tempImportance <- tempImportance[newOrder,] # orders complete
        # importance output considering previous arrangement.
        
        retention <- floor(nrow(tempImportance) * innerRetentionProportion)
        # Number of features retained as a result of eliminating a proportion
        # of the least important ones.
        cat(paste("Retained features (iteration: ",
                  as.character(iterations),
                  "; retention loop: ",
                  as.character(retentionLoop),
                  "): ",
                  as.character(retention),
                  ".\n",
                  sep = "")) # Message informing about the number of outer loop
        # iterations or number of times the entire data set has been used from
        # its entirety, the number of times features have been eliminated
        # starting from all the features until there is only one feature left
        # in the data set and the number of features that are currently left in
        # the data set.
        
        featsRetained <- rownames(tempImportance)[1:retention] # Actual name
        # of all features that will be kept for the next loop eliminating
        # iteration.
        
        featsEliminated <- rownames(tempImportance)[!(rownames(tempImportance) %in% featsRetained)]
        # All features that have been eliminated during this iteration.
        
        allFeatsEliminated <- append(x = allFeatsEliminated,
                                     values = featsEliminated)
        # Eliminated features are added to the vector which keeps their names
        # and increases its length after each round that passes.
        
        # This previous vector is used to incorporate information derived from
        # it to the matrix shell that saves the information about features
        # importance. Features already eliminated will receive a value of zero
        # for each round they do not appear among the most important features
        # left, while these ones will receive a value of 1. These values will
        # be added to the previously mentioned matrix, which arranges the
        # compilation of zeros and ones per columns, which are named after each
        # feature.
        
        for(i in 1:length(featsRetained)){ # This loop iterates over the name
          # of all features that have not been deleted from the data matrix
          # due to its relatively higher importance.
          
          ifeatRetained <- featsRetained[i] # selects a specific feature.
          iindex <- which(is.na(hitArr[,ifeatRetained]))[1] # retrieves the
          # first empty column position attributed to that feature due to
          # the presence of NA values.
          hitArr[iindex, ifeatRetained] <- 1 # fills it with a 1.
        
        }
        
        for(i in 1:length(allFeatsEliminated)){ # This loop iterates over the
          # name of all features that have already been deleted from the data
          # matrix due to its relatively lower importance.
          
          ifeatEliminated <- allFeatsEliminated[i] # selects a specific feature.
          iindex <- which(is.na(hitArr[,ifeatEliminated]))[1] # retrieves the
          # first empty column position attributed to that feature due to
          # the presence of NA values.
          hitArr[iindex, ifeatEliminated] <- 0 # fills it with a 0.
        
        }
        
        tempDataSet <- tempDataSet[, featsRetained] # The data set is reduced
        # to those features that have been considered proportionally more
        # important.
        
        retentionLoop <- retentionLoop + 1 # The number of loops during which
        # feature elimination is produced is increased in one unit. This
        # number will determine when inner loop iterations shall stop.
        
      }
      
      iterations <- iterations + 1 # If there is one feature left in the
      # data set, the number of outer loop iterations is increased by one.
      # This number will determine when outer loop iterations shall stop.
      
    }
    
    # After all importance measurement calculation and importance summation has
    # been performed, it is time to complete the importance hit matrix with
    # ones and zeros so a binomial test can be computed for all features.
    # The number of events per feature has to be equal and no NA values should
    # appear for this computation to be run.
    
    # To do so, NA values associated to each feature will be replaced by
    # a proportioned number of hits and zeros doing resampling from the actual
    # values accumulated for each feature.
    
    for(i in 1:ncol(hitArr)){ # For each feature represented in the importance
      # hit matrix.
      
      ibool <- is.na(hitArr[,i]) # A boolean vector that corresponds to those
      # column positions for which there are no attributed values.
      ihit <- hitArr[!ibool, i] # A boolean vector that corresponds to those
      # column positions for which there are attributed values (either 0 or 1).
      isample <- sample(x = ihit, size = sum(ibool), replace = TRUE)
      # executes resampling over present feature values as many times as
      # the number of NA values associated to it so they can be replaced.
      hitArr[ibool, i]<- isample # replaces NA values using the resampling
      # result.
    }
    
    # To compute the binomial test and be able to obtain a p-value from a set
    # of zeros and ones, the probability of a hit has to be calculated. It can
    # be computed as a proportion from the total number of events that are
    # being considered by the array.
    
    hitArrSize <- prod(dim(hitArr)) # Size of the array.
    totalHits <- sum(hitArr, na.rm = TRUE) # Sum of all hits.
    pHit <- totalHits / hitArrSize # Probability of an event being a hit.
    
    # Only when the number of hits is significantly high the obtained p-value
    # will be significant too, i. e., if the number of hits is significantly
    # low, the statistical test will not provide a significant p-value.
    
    allPValues <- apply(X = hitArr, # Importance hit array.
                        MARGIN = 2, # Array calculations will be executed by
                        # columns.
                        FUN = function(x) binom.test(x = sum(x, na.rm = TRUE),
                                                     # Number of hits per feature.
                                                     n = length(x),
                                                     # Total number of events per feature.
                                                     p = pHit, # Hit probability.
                                                     alternative = "greater")$p.value
                                                     # Only a high number of
                                                     # hits will contribute to
                                                     # greater significance.
                        ) # computes a p-value from a binomial test for each
    # feature based on random forest importance results.
    
    importanceMatrix[, "binomial.pvalue"] <- allPValues[rownames(importanceMatrix)]
    # incorporates binomial p-values to importance results matrix.
    importanceMatrix[, "adj.pvalue"] <- p.adjust(p = allPValues[rownames(importanceMatrix)],
                                                 method = pAdjMethod)
    # adjusts p-values according to scripted adjustment method. 
    
    orderMatrix <- order(importanceMatrix[,"importance"], decreasing = TRUE)
    # computes an order of the importance results matrix rows by decreasing raw
    # importance value. Each row represents a feature.
    importanceMatrix <- importanceMatrix[orderMatrix,]
    # orders matrix rows according to previous computation.
    
    boolOrdGreat <- importanceMatrix[,"importance"] >= 0 # retrieves those
    # indices representing the features with a positive raw importance from
    # the importance overall results matrix.
    boolOrdLess <- importanceMatrix[,"importance"] < 0 # retrieves those
    # indices representing the features with a negative raw importance from
    # the importance overall results matrix.
    
    sumGreat <- sum(importanceMatrix[boolOrdGreat,"importance"])
    # sums all positive raw importance values.
    sumLess <- sum(importanceMatrix[boolOrdLess,"importance"])
    # sums all negative raw importance values.
    
    importanceMatrix[boolOrdGreat,"normalized.importance"] <-
      importanceMatrix[boolOrdGreat,"importance"] / sumGreat # calculates
    # normalized importance values considering those positively important
    # features only.
    importanceMatrix[boolOrdLess,"normalized.importance"] <-
      -(importanceMatrix[boolOrdLess,"importance"] / sumLess) # calculates
    # normalized importance values considering those negatively important
    # features only.
    
    importanceMatrix[,"standard.cumulative"] <- cumsum(importanceMatrix[,"normalized.importance"])
    
    # A feature set based on their relatively higher importance can be computed
    # by using the elbow method on positive importance values:
    
    positiveImp <- rownames(importanceMatrix)[importanceMatrix[,"importance"] >= 0]
    # retrieves feature names that have a raw positive importance value.
    nPositiveImp <- length(positiveImp) # calculates number of features with a
    # raw positive importance value.
    slopeCalculation_x <- c(1, nPositiveImp) # A vector of length 2. Its first
    # element equals 1, and its second element equals the number of features
    # with positive importance. These numbers represent the indices of the most
    # important feature and the least important feature with positive
    # importance. They will be used to calculate a line which will be utilized
    # to computationally apply the elbow method.
    slopeCalculation_y <- c(max(importanceMatrix[positiveImp, "importance"]),
                            min(importanceMatrix[positiveImp, "importance"]))
    # retrieves the greatest importance value and the least importance value
    # considering those features with a positive importance result only.
    
    slope <- (slopeCalculation_y[2] - slopeCalculation_y[1]) /
      (slopeCalculation_x[2] - slopeCalculation_x[1]) # computes a rate of
    # differences from the values pairs of the previous vectors.
    
    interceptorsPositive <- c() # This vector will hold the interceptor value
    # after making a line with the previous calculated slope as a differences
    # rate through all the points that each feature with a positive importance
    # represents.
    
    # Now, all positive importance features will be represented in a
    # two-dimensional plot. The horizontal axis will just describe features
    # indices after ordering all these features from the most important one
    # to the least important one. The vertical axis will represent their
    # importance measurement.
    
    # If a line with the slope that has been calculated before is forced to
    # pass through each point, an interceptor can be calculated for each
    # feature. There where the interceptor minimizes can be considered the
    # "elbow" from the elbow method. That feature along with all those having a
    # higher importance value can be grouped as a set of heuristically
    # important features.
    
    for(j in 1:nPositiveImp){ # Loop that iterates over all indices of
      # features with a positive importance value.
      
      jx <- j # Feature index.
      jy <- importanceMatrix[positiveImp[j], # Name of the currently selected
                             # feature based on its index.
                             "importance"] # retrieves importance value of the
      # selected feature.
      jinterceptor <- jy - jx * slope # calculates interceptor associated to a
      # point representing a feature from a constant line slope.
      interceptorsPositive <- append(x = interceptorsPositive, values = jinterceptor)
      # adds interceptor result to interceptors vector.
      
    }
    
    importanceMatrix[,"y0.elbow.method"] <- c(interceptorsPositive,
                                              rep(NA, (nrow(importanceMatrix) - nPositiveImp)))
    # incorporates interceptors values to importance results matrix.
    importanceMatrix[,"is.elbow.method"] <- c(as.numeric(1:nPositiveImp < 
                                                           which.min(interceptorsPositive)),
                                              rep(0, (nrow(importanceMatrix) - nPositiveImp)))
    # computes boolean vector which indicates if each feature can be considered
    # into the group of heuristically important features based on the elbow
    # method on importance values.
    
    if(innerIsElbowPlot){ # If it was indicated that importance results and
      # elbow method results should be portrayed as a plot, this part of the
      # code will be executed.
      
      if(length(unique(as.character(innerGroups))) == 2){ # If the number of
        # classes in classification analysis equals 2, this portion of code is
        # run.
        
        # This entire function has the possibility of being used to obtain
        # importance results by pairs, contrasting each class to the rest of
        # them if the number of classes is greater than 2, or just gives the
        # possibility of obtaining importance results by comparing all classes
        # at once in random forests.
        
        # The following lines were written to address importance comparisons
        # in pairs. The class named "Other" clusters all samples from the
        # different classes that do not coincide with the one that is being
        # compared to the rest combined. "Other" represents the tag conceded
        # by this function to all samples from the different classes to
        # compute paired comparisons.
        
        # In a looped fashion, there will be as many importance plots as
        # paired comparisons if they have been requested.
        
        if(unique(as.character(innerGroups))[1] == "Other" &
           unique(as.character(innerGroups))[2] != "Other"){ # For paired
          # comparisons, from a factor object, "Other" class tag can be the
          # first one or the last one after alphabetical default ordering.
          # That is why this conditional and the next one try to detect if
          # any of those two tags, regardless their order, correspond to
          # that given string. In other words, if paired comparisons have
          # been requested, any of these two conditionals will be run.
          
          # The reason why this code layout has been chosen is because it
          # guarantees that the "Other" class tag will correspond to
          # "secondClass" variable content, although other code expressions
          # could have been used.
          firstClass <- unique(as.character(innerGroups))[2]
          secondClass <- unique(as.character(innerGroups))[1]
        }
        if(unique(as.character(innerGroups))[2] == "Other" &
           unique(as.character(innerGroups))[1] != "Other"){
          firstClass <- unique(as.character(innerGroups))[1]
          secondClass <- unique(as.character(innerGroups))[2]
        }
        # Otherwise, if results do not derive from multiple paired importance
        # comparisons but there are only two classes, this code chunk will be
        # run:
        if(unique(as.character(innerGroups))[1] != "Other" &
           unique(as.character(innerGroups))[2] != "Other"){
          firstClass <- unique(as.character(innerGroups))[1]
          secondClass <- unique(as.character(innerGroups))[2]
        }
        
        elbowPlotTitle <- paste("Elbow method on feature importance ('",
                                firstClass,
                                "' vs. '",
                                secondClass,
                                "')",
                                sep = "") # composes plot title from two class
        # tags.
        
      }else{ # In multiclass classification problems, this next string will be
        # used as elbow method results plot title:
        
        elbowPlotTitle <- paste("Elbow method on feature importance",
                                sep = "")
        
      }
      
      plot(x = 1:sum(importanceMatrix[,"is.elbow.method"]),
           # for the horizontal axis, a set of indices representing those
           # features with a positive importance value is chosen as plot
           # data.
           y = importanceMatrix[1:sum(importanceMatrix[,"is.elbow.method"]),
                                "normalized.importance"],
           # Corresponding normalized importance value for the selected
           # features is used as the vertical axis values.
           xlim = range(pretty(c(1,nrow(importanceMatrix)))),
           # establishes horizontal axis limits. Right limit considers the
           # number of all features in the analysis.
           xaxt = "n", # No axis tags marks or values will be plotted.
           ylim = range(pretty(importanceMatrix[,"normalized.importance"])),
           # establishes vertical axis limits based on normalized importance
           # values range.
           xlab = "Feature index", # names horizontal axis tag.
           ylab = "Random forest importance value", # names vertical axis tag.
           main = elbowPlotTitle, # gives plot title.
           type = "l", # Type of plot: continuous line.
           col = "red", # Color of the line.
           lwd = 2.5 # Factor that alters line width.
           )
      
      axis(side = 1, # Side on which a series of tags is wanted to be plotted:
           # conventional horizontal axis.
           at = unique(as.integer(pretty(c(1,nrow(importanceMatrix))))),
           # Positions on which tags are wanted to appear.
           labels = unique(as.integer(pretty(c(1,nrow(importanceMatrix)))))
           # Axis labels. This vector has to be equal in length to previous
           # argument vector.
           )
      
      points(x = sum(importanceMatrix[,"is.elbow.method"]):nrow(importanceMatrix),
             # Indices representing those features that do not have enough
             # importance based on the elbow method.
             y = importanceMatrix[sum(importanceMatrix[,"is.elbow.method"]):nrow(importanceMatrix),
                                  "normalized.importance"],
             # Importance value of those features without enough importance to
             # be selected by the elbow method.
             type = "l", # Type of points to be plotted: Joined by a line.
             col = "royalblue3", # Color of the line.
             lwd = 2.5 # Factor that determines line width.
             ) # Set of points to be added to previous plot.
      
      legend(x = "topright", # Wanted quadrant in which this legend will be
             # plotted.
             legend = c("Features selected",
                        "Features excluded"), # Legend text.
             inset = c(0.02, 0.04), # Factors used to move legend around
             # considering its quadrant default position.
             pch = 21, # Type of symbol wanted for the legend: They all will be
             # colored round points.
             col = "black", # Color of marker symbol perimeter.
             pt.bg = c("red", "royalblue3"), # Colors chosen for legend symbols.
             cex = 1.3 # Factor that modulates legend size.
             ) # incorporates a legend to the plot.
      
    }
    
    return(list("impArr" = importanceMatrix, 
                "hitArr" = hitArr)
           ) # Function output: A list containing the importance overall matrix
    # results and the importance hit array as two different objects.
    
  }
  
  if(length(unique(as.character(groups))) == 2){ # If there are only two groups
    # in the classification problem, it cannot be possible to apply the
    # binary contrasts approach multiple times, so one of the main function
    # arguments is forced to equal "FALSE":
    isBinaryApproach <- FALSE
  }
  
  if(!(isBinaryApproach)){ # If multiple binary comparisons approach for
    # multiclass classification problems has not been selected, this part of
    # the code will be executed.
    
    cat("\nNon-binary feature selection and classification approach chosen.\n")
    
    multiclassImpMatrix <- rrfIterator(innerDataSet = dataset, # Data set.
                                       innerGroups = groups, # Vector  of
                                       # classes tags (as a factor object).
                                       innerTreeNumber = treeNumber, # Number
                                       # of decision trees per random forest.
                                       innerRetentionProportion = retentionProportion,
                                       # Proportion of features that will be
                                       # retained in each elimination round.
                                       innerRetentionLoops = retentionLoops,
                                       # Number of times the data set will be
                                       # processed from its entirety to
                                       # obtain importance measurements by
                                       # feature elimination.
                                       innerIsElbowPlot = isElbowPlot
                                       ) # obtains importance results for a
    # given data set based on the iterative application of random forest.
    
    if(chosenFeatsClassMethod == "elbow"){ # If the method to select a group of
      # features has been determined to be the elbow method, this part of
      # the code will be used. The resulting features will be used to design a
      # last random forest classifier.
      chosenFeatures <-
        rownames(multiclassImpMatrix[["impArr"]])[as.logical(multiclassImpMatrix[["impArr"]][,"is.elbow.method"])]
    }
    if(chosenFeatsClassMethod == "p"){ # This conditional will be executed if
      # the criterion considered to select a subset of the most important
      # features is their associated p-value after the results of this
      # algorithm.
      chosenFeatures <-
        rownames(multiclassImpMatrix[["impArr"]])[multiclassImpMatrix[["impArr"]][,"is.elbow.method"] <= chosenPValue]
    } # If it is just wanted to consider a number of the most important
    # features, then this conditional is executed. "kFeats" variable is needed
    # to indicate the number of important features that are wanted to be
    # retrieved.
    if(chosenFeatsClassMethod == "k"){
      chosenFeatures <- rownames(multiclassImpMatrix[["impArr"]])[1:kFeats]
    }
    
    finalRF <- randomForest(x = dataset[, chosenFeatures],
                            y = groups,
                            ntree = treeNumber,
                            importance = TRUE,
                            proximity = TRUE
                            ) # It generates a last random forest classifier
    # only with those features that have been selected as important. It should
    # optimize the classification procedure. It has been indicated that
    # importance and random forest sample proximity will be calculated and
    # output.
    
    return(list("importance.matrix" = multiclassImpMatrix,
                "selected.features" = chosenFeatures,
                "RF.classifier" = finalRF)
           ) # List of results having all importance results, a set of
    # important features considered by a selected criterion and the final
    # random forest classifier.
    
  }else{ # This part of the code will be run if it has been established that
    # binary contrasts will be executed as an approach for a multiclass problem.
    
    cat("\nBinary feature selection and classification approach chosen.\n")
    
    listOfBinaryClassifiers <- list() # creates an empty list to store
    # importance results for binary contrasts.
    allChosenFeatures <- c() # constitutes a vector that will be iteratively
    # filled with those features considered important to participate in
    # a correct classification between the different classes and the rest of
    # them combined under a unifying class tag. The set of important features
    # retrieved each time will take into account the criterion determined by
    # the argument "chosenFeatsClassMethod".
    
    for(i in 1:nlevels(groups)){ # Loop that iterates over the indices of
      # class tags. It will be used to compare each class to the rest of
      # samples belonging to different classes at the same time.
      
      iclass <- unique(as.character(groups))[i] # selects one class.
      igroups <- as.character(groups) # converts classes vector into a
      # character vector object.
      igroups[igroups != iclass] <- "Other" # concedes class tag "Other" to
      # all elements that do not belong to the previously chosen class.
      igroups <- as.factor(igroups) # transforms modified class vector into
      # a factor object.
      
      cat(paste("\nChosen class in binary procedure: '", iclass, "'.\n", sep = ""))
      
      ibinaryImpMatrix <- rrfIterator(innerDataSet = dataset, # Data set.
                                      innerGroups = igroups, # Classes tags
                                      # vector.
                                      innerTreeNumber = treeNumber, # Number of
                                      # trees per random forest.
                                      innerRetentionProportion = retentionProportion,
                                      # proportion of features that will be
                                      # kept after each feature elimination
                                      # round.
                                      innerRetentionLoops = retentionLoops,
                                      # Number of times the data set will be
                                      # used from its entirety to extract
                                      # information about feature importance
                                      # by iteratively deleting all features
                                      # except one and starting back again.
                                      innerIsElbowPlot = isElbowPlot
                                      # determines if an elbow method results
                                      # representation is desired for each
                                      # paired contrast.
                                      ) # generates importance results
      # associated to a binary contrast.
      
      if(chosenFeatsClassMethod == "elbow"){ # is executed if the elbow method
        # is the wanted criterion to retrieve a set of important features for
        # each one of the binary contrasts in a multiclass problem.
        ichosenFeatures <-
          rownames(ibinaryImpMatrix[["impArr"]])[as.logical(ibinaryImpMatrix[["impArr"]][,"is.elbow.method"])]
      }
      if(chosenFeatsClassMethod == "p"){ # is executed if the statistical
        # significance calculated from the binomial test generated data
        # is the wanted criterion to retrieve a set of important features for
        # each one of the binary contrasts in a multiclass problem.
        ichosenFeatures <-
          rownames(ibinaryImpMatrix[["impArr"]])[ibinaryImpMatrix[["impArr"]][,"is.elbow.method"] <= chosenPValue]
      }
      if(chosenFeatsClassMethod == "k"){ # is executed if just a set of the
        # most important features by a decreasing importance order is the
        # wanted criterion to retrieve a set of important features for each one
        # of the binary contrasts in a multiclass problem.
        ichosenFeatures <- rownames(ibinaryImpMatrix[["impArr"]])[1:kFeats]
      }
      
      allChosenFeatures <- append(x = allChosenFeatures, values = ichosenFeatures)
      # The previous retrieved set of features is added to a vector that will
      # gather all those features that had importance in classifying the
      # samples from different classes in binary contrasts.
      
      ibinaryRF <- randomForest(x = dataset[,ichosenFeatures],
                                # uses a data set only with the selected
                                # important features.
                                y = igroups, # binary tags classes vector.
                                ntree = treeNumber, # Number of decision trees
                                # per random forest.
                                importance = TRUE, # computes feature
                                # importance measurement.
                                proximity = TRUE # computes random forest
                                # sample similarity.
                                )
      
      listOfBinaryClassifiers[[iclass]] <- list("importance.matrix" = ibinaryImpMatrix,
                                                "chosen.features" = ichosenFeatures,
                                                "RF.classifier" = ibinaryRF
                                                ) # gathers all results from
      # the binary classification in a single list.
      
    }
    
    allChosenFeatures <- unique(allChosenFeatures) # deletes feature duplicates
    # from the character vector.
    
    finalRF <- randomForest(x = dataset[, allChosenFeatures],
                            y = groups,
                            ntree = treeNumber,
                            importance = TRUE,
                            proximity = TRUE
                            ) # executes a final classifier generation
    # considering all the features without repetition that have shown to be
    # important in any of the binary comparisons.
    
    return(list("binary.classifiers" = listOfBinaryClassifiers,
                "selected.features" = allChosenFeatures,
                "RF.classifier" = finalRF)
           ) # A list of results returning the importance matrix results of
    # all binary comparisons, what features have been used to create the
    # last classifier and the final random forest classifier.
    
  }
  
}

## 7) reducedDescription:

# This function creates a set of substrings from a character vector having
# protein descriptions as elements which information has been obtained from
# UniProt. Each substring consists of all characters that precede a series of
# UniProt abbreviations that contain information about the species, the protein
# associated gene, its sequence version, etc. the protein has to do with.

reducedDescription <- function(fullDescription){
  
  # The given input for this function consists of a matrix having in its second
  # column the character vector with the proteins descriptions.
  
  re <- regexpr(".* OS=", fullDescription[,2]) # retrieves all characters
  # indices that precede the first UniProt description abbreviation code.
  ss <- substr(x = fullDescription[,2],
               start = re,
               stop = attr(x = re, which = "match.length") - 4) # creates
  # a substring out of the original text by considering previous retrieved
  # character indices.
  best <- paste(ss, " (", rownames(fullDescription), ")", sep = "") # adds
  # UniProt protein code name between brackets at the end of the reduced
  # description.
  names(best) <- rownames(fullDescription) # names elements of a vector using
  # the row names from the given input matrix.
  return(best)
  
}

## 8) multiBox2:

# This function creates a boxplot representing the levels of abundance of
# a given feature for multiple classes or groups and incorporates a
# statistical significance relation between groups if given.

multiBox2 <- function(boxValues, # A series of numerical data.
                      featName, # Name of the analyzed feature.
                      groups, # Character vector representing sample groups
                      # corresponding to numerical data order.
                      pvals, # A series of p-values representing statistical
                      # significance between pairs of groups. Vector p-values
                      # order is specific and determines this function correct
                      # performance. It depends on alphabetical group name
                      # order and all possible paired group combinations.
                      onlySig = TRUE, # A boolean variable that specifies if
                      # the only kinds of relations between groups that are
                      # wanted to be plotted imply statistical significance.
                      chosenPVals = FALSE, # determines if there is only a
                      # group of features which statistical significance
                      # is wanted to be represented.
                      whatPVals, # If "chosenPVals = TRUE", using this
                      # argument what statistical significant comparisons are
                      # wanted to be plotted can be selected.
                      textCex = 1.5, # Factor that controls statistical
                      # significance text symbol size.
                      main2, # Plot title.
                      isNorm = TRUE # indicates if values have been previously
                      # normalized to adequately modify vertical axis label
                      # text.
                      ){
  
  boxDataF <- data.frame(boxValues, groups) # builds a data frame containing
  # the set of feature values in one column and the ordered sample class
  # belonging in a second column, with which data are orderly corresponded to
  # the specific classes.
  colnames(boxDataF) <- c(featName, "classes") # renames previous data frame
  # columns using the feature name and "classes" string.
  
  # This function is able to indicate using parallel lines to the horizontal
  # plot axis statistical significance relations between groups, using the
  # conventional significance asterisk symbols. To plot these symbols, line
  # coordinates are calculated based on the number of groups, groups position
  # on the plot and reference height based on plot elements height.
  
  rights <- 1:(length(unique(groups))) + 1/(length(unique(groups))*4)
  # creates a set of points for the plot horizontal axis. This set of points
  # corresponds to a group of consecutive integers starting from one, having
  # their consecutive elements distanced by one unit, and ending with a number
  # equal to the total number of classes, all of them having a constant
  # proportion added. This proportion equals one divided to a product, which
  # equals the number of classes times 4. This implies that this constant
  # proportion changes with the number of classes and is reduced with the
  # increase of total classes.
  rights[length(rights)] <- length(rights) - 1/(length(rights)*4)
  # The last previous numerical vector element is changed to the total number
  # of classes less than the previous defined constant proportion.
  lefts <- 1:(length(unique(groups))) - 1/(length(unique(groups))*4)
  # calculates a set of similar horizontal positions to the last calculated
  # vector positions but, in this case, displaced using the computed proportion
  # by subtracting it instead of by adding it to the group of integer indices.
  lefts[1] <- 1 + 1/(length(lefts)*4) # corrects first vector element by adding
  # the constant calculated proportion instead of subtracting it.
  
  firstPretty <- range(pretty(boxValues))[2] # ranges a set of values from a
  # minimum to a maximum based on input feature values and retrieves only the
  # maximum value.
  secondPretty <- range(pretty(boxValues))[1] # ranges a set of values from a
  # minimum to a maximum based on input feature values and retrieves only the
  # minimum value.
  difPretty <- firstPretty - secondPretty # calculates difference between
  # maximum and minimum retrieved values that have just been retrieved.
  
  # The next loop will create all possible paired combinations of class indices
  # and class tags without repetition. In the case of indices combinations, the
  # first index will always be the smallest number. This will help retrieving
  # the combination of positions that have to be taken into account to plot
  # the lines that compare groups on the plot and establish a relation of
  # statistical significance.
  
  indComb <- list() # A list that will compile paired class tags indices
  # combinations without repetition.
  namesIndComb <- c() # A vector that will compile paired class tags names
  # combinations.
  k <- 1 # An index that will be iteratively updated to order the paired
  # indices combinations list.
  for(i in 1:(length(unique(groups)) - 1)){ # Outer loop that iterates over
    # class tags indices except the last one in increasing order.
    for(j in (i + 1):length(unique(groups))){ # Inner loop that iterates over
      # class tags indices except the first one in increasing order.
      indComb[[k]] <- c(i, j) # creates a pair of class tags indices and adds
      # it to a list. To index it to the list and sequentially add the created
      # pair of indices, "k" index is used.
      tempPName <- paste(unique(groups)[i], unique(groups)[j], sep = "-")
      # creates paired class tags string based on using outer and inner
      # loops indices.
      namesIndComb <- append(x = namesIndComb, values = tempPName) # adds
      # paired class tag combination to the previously created vector.
      k <- k + 1 # "k" index is updated at the end of the inner loop.
    }
  }
  names(indComb) <- namesIndComb # names combinations indices list with
  # the character vector containing paired class tags combinations.
  
  sigHeights <- c() # A vector which purpose is to compile the plot heights of
  # the horizontal lines that establish statistical significance relations
  # between the boxplot depicted groups.
  textHeights <- c() # A vector compiling the text that will appear on the
  # previous lines to show the kind of statistical significance power between
  # the two related groups considering conventionally used significance symbols
  # ("n.s.", "*" if p-value < 0.05, "**" if p-value < 0.01, "***" if
  # p-value < 0.001, or "****" if p-value < 0.0001).
  
  if(onlySig){ # If it has been determined that only significant contrasts
    # will be indicated through "onlySig" argument when running this
    # user-defined function, this conditional will be executed.
    whichPVals <- which(pvals <= 0.05) # creates a boolean vector that
    # indicates which p-values are less than 0.05 (which represents the
    # significance threshold). 
  }
  if(chosenPVals){ # If it has been determined that only a set of specific
    # contrasts will be marked as significant or not through the "chosenPVals"
    # argument, this conditional will be executed.
    whichPVals <- whatPVals
  }
  
  p <- 1 # "p" defines a variable to compute multiple plot heights destined
  # to represent statistical significance comparison lines between groups.
  q <- 2 # "q" defines a variable that will generate a set of values slightly
  # greater than those associated to "p" to be used to plot text over the
  # previously mentioned statistical significance contrasts lines between
  # groups. "q" will be used to plot the kind of significance degree between
  # groups. This is why "p" is slightly less than "q", so the text can appear
  # in a subtly higher position than the contrast line it refers to.
  for(i in 1:length(whichPVals)){ # A loop that iterates over all p-values
    # indices from the p-values vector.
    newSig <- max(boxValues) + (p*(difPretty/50)) # sums maximum value from
    # input data and a product together. The product equals the variable "p"
    # times the result of a rate. This rate equals the difference between the
    # maximum value and the minimum value of a range calculated before
    # considering all input data divided by a constant which equals 50.
    # This value increases as "p" does. It will represent a statistical
    # significance line height that will be used to plot it. For each
    # significant contrast between groups, a line height value will be produced. 
    sigHeights <- append(x = sigHeights, values = newSig) # adds recently
    # computed sum to previously created vector.
    newText <- max(boxValues) + (q*(difPretty/50)) # This calculation is
    # equivalent to the kind of computation that is being executed to obtain
    # a height value associated to a statistical significance contrast line.
    # In this case, the variable "q" is being used instead of "p", which is
    # greater than the latter by one unit. This implies that this sum will
    # be slightly greater than the one considering "p". This calculated height
    # value will be used to plot a text over the line which will refer to
    # the statistical significance degree for this comparison.
    textHeights <- append(x = textHeights, values = newText) # adds height value
    # to numerical vector.
    p <- p + 2 # "p" variable is increased by 2 units.
    q <- q + 2 # "q" variable is increased by 2 units.
  }
  
  finalPretty <- max(textHeights) # stores maximal text height in a separate
  # variable.
  
  boxDataF$classes <- factor(boxDataF$classes , levels = unique(groups))
  # converts data frame column having the samples classes into a factor object.
  # Plus, the level order is determined through the variable that stores
  # samples class tags. After computing what class tags without repetition
  # are and using them as the factor levels, boxplot classes will be depicted
  # in the same order. Therefore, sample classes ordering in the input data
  # array or data frame affects the corresponding classes tags vector and,
  # subsequently, boxplot classes depiction disposition.
  
  if(isNorm){ # After indicating with a boolean variable if the used values
    # from the data array are normalized or not, the corresponding part of
    # this conditional will be run, which will affect the vertical axis label
    # text content by indicating if feature values have been normalized or not.
    boxplotYLab <- paste(featName, " normalized levels", sep = "")
    # stores string to represent that data values have been previously
    # normalized.
  }else{
    boxplotYLab <- paste(featName, " raw levels", sep = "")
    # stores string to represent that data values have not been previously
    # normalized.
  }
  
  boxplot(boxDataF[,1] ~ boxDataF[,2], # uses formula format to make the
          # boxplot.
          data = boxDataF, # Data frame from which values and class tags are
          # obtained.
          pch = 21, # Factor that affects boxplot point size.
          bg = "grey", # colors boxplot points.
          xaxt = "n", # impedes any horizontal axis numerical labels to be
          # plotted.
          xlab = "Groups", # names horizontal axis.
          ylab = boxplotYLab, # names vertical axis.
          main = main2, # Plot title.
          ylim = range(pretty(c(min(boxValues), finalPretty))) # establishes
          # vertical axis upper and lower plot limits.
          ) # generates boxplot image.
  
  axis(side = 1, at = 1:length(unique(groups)), labels = unique(groups))
  # plots class tags under each corresponding group box image.
  
  if(length(whichPVals) >= 1){ # If the boolean vector indicating which
    # p-values from their numerical vector have to be represented as
    # statistical significance symbols has a positive length, this conditional
    # will be executed.
    for(i in 1:length(whichPVals)){ # This loop will iterate over all p-values
      # boolean vector elements indices.
      
      if(length(unique(groups)) == 2){ # If there are only 2 different class
        # tags in the analysis, this code chunk will be run.
        tempStart <- 1 # Variable "tempStart" will constantly equal 1. This
        # variable refers to the class box position from left to right that
        # has to be considered to be the starting point of statistical
        # significance horizontal contrast lines. Since there are 2 groups,
        # there can only be one comparison or contrast between groups, so
        # there can be only one kind of line matching the two groups.
        tempEnd <- 2 # Variable "tempEnd" will constantly equal 2. This
        # variable refers to the class box position from left to right that
        # has to be consider to be the ending point of statistical
        # significance horizontal contrast lines.
      }
      ti <- whichPVals[i] # iteratively stores a boolean value from the p-values
      # boolean vector.
      if(length(unique(groups)) > 2){ # If the number of classes is greater
        # than 2, this conditional will be executed.
        
        # It is assumed that p-values come from a vector named using a specific
        # pattern that comes from another user-defined function created to
        # calculate statistical significance p-values between multiple groups
        # through an ANOVA post hoc-test. Name tags have an established
        # character pattern and the result from this function are derived from
        # them.
        
        # Next two code commands make sure two possible substrings are not
        # present anymore in the given p-values vector elements names.
        
        newnames <- gsub(pattern = "adj.pval.", # pattern to be found in the
                         # given string.
                         replacement = "", # Pattern that will substitute the
                         # last pattern indicated.
                         x = names(ti), # String from which a pattern will be
                         # found and replaced by another one.
                         fixed = TRUE # Argument used to determine searched
                         # pattern in the given string is exactly it and
                         # it does not refer to a regular expression.
                         ) # substitutes a substring by another one.
        newnames <- gsub(pattern = "pval.", # Pattern to be found in the
                         # given string.
                         replacement = "", # Pattern that will substitute the
                         # last pattern indicated.
                         x = newnames, # String from which a pattern will be
                         # found and replaced by another one.
                         fixed = TRUE # Argument used to determine searched
                         # pattern in the given string is exactly it and
                         # it does not refer to a regular expression.
                         ) # substitutes a substring by another one.
        
        group1 <- strsplit(x = newnames, # Input character object.
                           split = "-", # indicates character substring from
                           # which to induce string splitting.
                           fixed = TRUE)[[1]][1]
        # splits character vector elements or a single character object into
        # a character vector. It splits them considering a character or set
        # of them from which characters elements will be divided into multiple
        # strings. The output object is a list, where each element is a
        # character vector.
        
        # Considering the known origin of the input data that are used into
        # this function, it is known that from this string splitting there
        # will be two substrings produced, each one being a class tag. The
        # next code line stores the second class tag.
        
        group2 <- strsplit(x = newnames,
                           split = "-",
                           fixed = TRUE)[[1]][2]
        # generates class tag out of string splitting from p-values numerical
        # named vector labels.
        
        orderMe <- c() # creates empty vector.
        orderMe <- append(x = orderMe, # indicates vector in which values are
                          # wanted to be added to.
                          values = which(unique(groups) == group1) # calculates
                          # index from the class tags vector without any class
                          # tag repetition that corresponds to the first
                          # retrieved group out of the string splitting
                          # computation. This index corresponds to a box
                          # position from the boxplot image.
                          ) # adds class tag index out of class tags vector
        # into numerical compiling vector.
        orderMe <- append(x = orderMe, # indicates vector in which values are
                          # wanted to be added to.
                          values = which(unique(groups) == group2) # calculates
                          # index from the class tags vector without any class
                          # tag repetition that corresponds to the second
                          # retrieved group out of the string splitting
                          # computation. This index corresponds to a box
                          # position from the boxplot image.
                          ) # adds class tag index out of class tags vector
        # into numerical compiling vector.
        orderMe <- orderMe[order(orderMe)] # orders class tags indices vector
        # considering numerical increasing order. It just stores 2 indices
        # which will now be numerically ordered.
        
        tempStart <- orderMe[1] # retrieves first class index from previous
        # vector. It indicates which is the class box position to be considered
        # to start plotting a statistical significance contrast line.
        tempEnd <- orderMe[2] # retrieves first class index from previous
        # vector. It indicates which is the class box position to be considered
        # to end plotting a statistical significance contrast line.
      }else{ # If the total number of classes in the analysis is 2, the
        # previous vector created to compile class tags indices is now
        # generated to have the only 2 class tags in alphabetical order:
        orderMe <- unique(groups)[order(unique(groups))]
      }
      
      tempStart <- rights[tempStart] # overwrites the first class tag index
      # variable to store the original specific horizontal axis position on
      # the plot to create the line representing a statistical significance
      # comparison between the two chosen groups.
      tempEnd <- lefts[tempEnd] # overwrites the second class tag index
      # variable to store the ending specific horizontal axis position on the
      # plot to create the line representing a statistical significance
      # comparison between the two chosen groups.
      tempMean <- (tempStart + tempEnd)/2 # calculates mean horizontal axis
      # position from the two line vertices retrieved before.
      littleSum <- difPretty/100 # Hundredth part of a distance obtained from
      # the vertical axis values range difference of the plot. This distance
      # is then related to a height obtained out of the vertical axis values.
      segments(x0 = tempStart,
               x1 = tempEnd,
               y0 = sigHeights[i]
               ) # draws a horizontal line between the two compared groups
      # considering the boxes on the plot representing them. This line is
      # parallel to the horizontal axis.
      segments(x0 = tempStart,
               y0 = sigHeights[i],
               y1 = sigHeights[i] - littleSum
               ) # draws a little line from the previously horizontal drawn
      # line to a subtly lower point, being completely parallel to the vertical
      # axis. This line starts at the rightest point of the horizontal
      # contrast line.
      segments(x0 = tempEnd,
               y0 = sigHeights[i],
               y1 = sigHeights[i] - littleSum
               ) # draws a little line from the previously horizontal drawn
      # line to a subtly lower point, being completely parallel to the vertical
      # axis. This line starts at the leftmost point of the horizontal
      # contrast line.
      tempPval <- pvals[ti] # orderly retrieves p-value from the numerical
      # p-values vector.
      
      # Next conditionals are ordered through code considering a decreasing
      # statistical significance value, to update the corresponding statistical
      # significance symbol until a certain point is reached, which depends on
      # the p-value that is being analyzed in each loop iteration. Depending on
      # the statistical significance degree of the analyzed classes contrast,
      # a symbol will be chosen and used as text to be represented on the plot.
      
      if(tempPval > 0.05){ # If the p-value is greater than the usual
        # statistical significance threshold 0.05, this conditional will be
        # run:
        tempText <- "ns"
      }
      if(tempPval <= 0.05){ # If the p-value is less than the usual
        # statistical significance threshold (0.05), this conditional will be
        # run:
        tempText <- "*"
      }
      if(tempPval <= 0.01){ # If the p-value is less than 0.01, this
        # conditional will be run:
        tempText <- "**"
      }
      if(tempPval <= 0.001){ # If the p-value is less than 0.001, this
        # conditional will be run:
        tempText <- "***"
      }
      if(tempPval <= 0.0001){ # If the p-value is less than 0.0001, this
        # conditional will be run:
        tempText <- "****"
      }
      
      text(x = tempMean, # Horizontal axis position on which the wanted text
           # will be plotted.
           y = textHeights[i], # Vertical axis position on which the wanted
           # text will be plotted.
           labels = tempText, # Text.
           cex = textCex # Factor which modulates text size.
           ) # plots statistical significance symbol text on the horizontal
      # contrast line between groups.
    }
  }
  
}

## 9) photoBox2:

# It makes use of the user-defined function "multiBox2" to plot a set of
# boxplots images as JPEG files where indicated.

photoBox2 <- function(is.already = TRUE, # Argument that takes a boolean
                      # value to indicate if the directory or folder where
                      # all images are going to be generated has been already
                      # created or not. If it is indicated that there is no
                      # folder to hold the images yet, this function will
                      # create it.
                      nameDir, # Name of the directory which will store all
                      # the produced JPEG images.
                      where, # Path to the folder that will hold all the
                      # images.
                      dat, # Data set containing all the values of features
                      # arranged by rows and the columns representing samples
                      # from which boxplots data will be obtained.
                      chosenQuality = 100, # Percentage that controls image
                      # quality. The greater it is, the lesser degradation of
                      # the image, although its compression will be lower.
                      chosenWidth = 600, # Factor that controls image width.
                      chosenHeight = 480, # Factor that controls image height.
                      pvalues, # P-values matrix where features are represented
                      # across rows and statistical contrasts significance is
                      # represented across columns. It is worth mentioning that
                      # these columns are ordered in a specific way considering
                      # classes tags alphabetical order and number of
                      # statistical comparisons considering the total number
                      # of groups. This ultimately depends on the results of
                      # an user-defined function that was first scripted to
                      # perform an ANOVA or Kruskal-Wallis post-hoc test using
                      # a data set and a vector indicating samples classes as
                      # input.
                      groups2, # Complete character orderly indicating the
                      # class tags each samples belong to considering data
                      # set rows.
                      is.ANOVA = TRUE, # A boolean value that was determined
                      # to be used to indicate if p-values were a result of
                      # an ANOVA post-hoc test or any other statistical test
                      # that provides statistical significance results
                      # considering three or more groups. If it equals "FALSE",
                      # then it would be considered that they come from a
                      # t-test or any other test for two classes only.
                      main3, # Boxplots titles as a character vector.
                      isNorm2 = TRUE # Argument that indicates if data set
                      # values have been normalized or not. It just changes
                      # vertical axis label to say if the used data values
                      # are normalized or if they just constitute raw values.
                      ){
  
  if(!(is.already)){ # If the folder that will hold all boxplot images has not
    # been created, this conditional will be run to create a directory.
    where <- paste(where, nameDir, "/", sep = "") # A string representing a
    # path to the given directory.
    dir.create(path = where) # creates folder on the given path in the previous
    # variable.
  }
  
  for(i in 1:nrow(dat)){ # Using the data matrix from the function input,
    # this loop iterates over its rows to retrieve a data row each time, which
    # will represent the necessary data to plot one boxplot image at a time.
    fileName <- paste(rownames(dat)[i], ".jpeg", sep = "") # creates the file
    # name of a boxplot image. Its name is based on the row names content the
    # input data set provides.
    filePath <- paste(where, fileName, sep = "") # pastes file name to boxplots
    # images directory path.
    jpeg(file = filePath, # Complete file path.
         quality = chosenQuality, # Chosen image quality.
         width = chosenWidth, # Chosen image width.
         height = chosenHeight # Chosen image height.
         ) # enables JPEG image file creation.
    
    if(is.ANOVA){ # If it has been indicated that p-values results have been
      # generated by a statistical analysis adapted to three or more classes,
      # this code chunk will be run.
      
      currentPValues <- pvalues[i,] # A series of p-values obtained from input
      # p-values matrix representing statistical significance between
      # pairs of groups. Vector p-values order is specific and
      # determines this function correct performance. It depends on
      # alphabetical group name order and all possible paired groups
      # combinations.
      
    }else{ # If it has been indicated that p-values results have been
      # generated by a statistical analysis adapted to only two classes, this
      # code chunk will be run.
      
      currentPValues <- pvalues[i] # A single p-value attributed to one feature
      # This p-value represents the statistical relation between two groups.
      
    }
    
    multiBox2(boxValues = dat[rownames(dat)[i],], # Data retrieved from the
              # first row of the input data set which has been indexed using
              # its row name.
              featName = rownames(dat)[i], # Orderly retrieved data set row
              # name.
              groups = groups2, # Complete character orderly indicating the
              # class tags each sample belongs to considering data set rows.
              pvals = currentPValues, # Selected p-values.
              onlySig = TRUE, # A boolean variable that specifies if the only
              # kinds of relations between groups that are wanted to be
              # plotted imply statistical significance. It is constantly
              # sustained as "TRUE" when this function is used.
              chosenPVals = FALSE, # determines if there is only a group of
              # features which statistical significance is wanted to be
              # represented. It is constantly sustained as "FALSE" when this
              # function is used.
              main2 = main3[i], # Plot title dragged from a character vector
              # as one of its ordered elements.
              isNorm = isNorm2 # indicates using a boolean value if data set
              # values have been normalized or not.
              ) # plots boxplot image.
    
    dev.off() # deactivates current open image device to be modified and
    # eventually creates image file.
  }
  
}

## 10) RFClassEx:

# It generates a bar plot from a random forest classifier object. Based on the
# features that the classifier has considered and the importance measurement
# random forest provides, this function creates a plot where a total importance
# value associated to each feature is represented as a bar. Plus, this total
# importance value also describes the relative importance to classify samples
# according to their class belonging using a color code. It provides an easily
# interpretable way to understand the power of data features to correctly
# classify samples.

# This kind of plot can be complemented with other representations that inform
# about the data sense, relatively positive or negative, each class aim to
# in comparison with other classes, to provide a deeper meaning of how correct
# classification is occurring. Thanks to data feature meaning interpretation,
# more knowledge can be drawn from machine learning classifiers to comprehend
# how they are reaching conclusions.

RFClassEx <- function(rfObject, # Random forest object from "randomForest"
                      # library usage.
                      trueClassesOrder, # A character vector to specify classes
                      # tags without repetitions that will establish which
                      # order classes occupy in relation to the color code used
                      # for bars and legend text.
                      classesColors, # A character vector indicating the
                      # ordered associated colors to their corresponding
                      # classes.
                      cexMain = 1.4, # Factor that modulates plot title size.
                      cexAxisNumbers = 1.3, # Factor that alters vertical axis
                      # numbers.
                      cexPlotAxis = 1.2, # Factor used to change axes names
                      # label size. 
                      tiltAngle = 40, # Number of degrees horizontal axis
                      # labels are tilted with respect to horizontal axis.
                      cexLabels = 1.2, # Number that modulates horizontal
                      # labels size.
                      cexLeg = 1.3, # Factor that controls plot legend size.
                      legInset = c(0.02, 0.02), # Numerical vector that
                      # displaces legend using two factors referred to the
                      # horizontal axis or vertical axis of the image.
                      adjText = c(1, 1) # A numerical vector that slightly
                      # displaces horizontal axis labels up, down, to the left
                      # or right.
                      ){
  
  if(!("package:randomForest" %in% base::search())){ # Conditional used to load
    # random forest library if it is not already loaded into R environment.
    catMe <- paste("randomForest library wasn't attached. Attaching library...\n")
    cat(catMe)
    library(randomForest) # loads R library.
    catMe <- paste("Done.\n")
    cat(catMe)
  }
  
  impArr <- rfObject[["importance"]][, trueClassesOrder] # retrieves importance
  # output matrix from the function input random forest object and establishes
  # a specific column order based on the input classes tags character vector,
  # used to index importance matrix columns.
  nFeatures <- nrow(impArr) # stores number of rows in the importance matrix,
  # which equals the number of features the classifier has considered while
  # being trained.
  
  barplotData <- t(impArr[order(colSums(t(impArr)), decreasing = TRUE),])
  # transposes importance matrix which originally had its columns representing
  # each class and its rows representing each feature, so its columns will now
  # be the features and its rows will be the classes. After that, a summation
  # across each column elements is computed, which will correspond to the
  # importance sum for a feature derived from an importance value for each
  # class. This set of results is decreasingly ordered, its ordering indices
  # are obtained from it and they are used to index the rows of the transposed
  # importance matrix. The final transposed importance matrix after indexing is
  # stored in a variable.
  
  percentTableValues <- t(apply(X = t(barplotData),
                                MARGIN = 1,
                                FUN = function(x) (x / sum(x)) * 100)
                          ) # As part of the output, a matrix will be produced
  # representing a proportion of what each specific value is equivalent to
  # after calculating this proportion from the summation of all classes
  # associated values. Therefore, a comparison of the importance value
  # between classes for each feature can be obtained.
  
  barplot(height = barplotData, # plots feature importance data in decreasing
          # order of importance.
          xaxt = "n", # avoids generating the horizontal axis numerical marks.
          col = classesColors, # character Vector of ordered color tags.
          ylab = "Class mean decrease accuracy (importance value)",
          # Vertical axis name label.
          ylim = range(pretty(c(0, max(t(apply(X = impArr,
                                               MARGIN = 1,
                                               FUN = cumsum)))))),
          # generates a numerical vector that establishes a range of a minimal
          # value and a maximal one to limit what is plotted as vertical axis,
          # guaranteeing an aesthetic plot image.
          yaxt = "n", # avoids generating the vertical axis numerical marks.
          main = "Class feature importance barplot", # Plot title.
          cex.main = cexMain, # Factor that controls plot title size.
          cex.lab = cexPlotAxis # Factor that controls axes labels tags.
          )
  axis(side = 2, # This value refers to vertical plot axis.
       at = pretty(c(0, max(t(apply(X = impArr, MARGIN = 1, FUN = cumsum))))),
       # generates a set of values to be represented as the numerical series
       # marks for the vertical axis.
       labels = pretty(c(0, max(t(apply(X = impArr, MARGIN = 1, FUN = cumsum))))),
       # generates a set of values to be represented as the numerical series
       # labels for the vertical axis.
       las = 1, # Factor that modifies text orientation. This values means that
       # text writing will be parallel to the horizontal axis.
       cex.axis = cexAxisNumbers # Factor that changes plotted text size.
       )
  
  spaceBetweenBars <- 0.2 # A value used in bar plot representation. It
  # determines the constant space between consecutive bars.
  barWidth <- 1 # A value used in bar plot representation. It determines the
  # the constant space bar width.
  root <- (spaceBetweenBars + barWidth/2) # Summation of the bar width
  # distance halved plus the distance between consecutive bars. It corresponds
  # to the middle horizontal point of the first bar and the horizontal position
  # the first feature label should occupy on the plot. From it (that is why
  # this variable is called "root"), the next middle points of the bars can be
  # deduced.
  centersDistance <- (spaceBetweenBars + barWidth) # Value that has to be
  # recursively summed to the previous initial value to obtain all tags
  # horizontal coordinates.
  xAxis <- c( root, (root + (1:(ncol(barplotData) - 1) * centersDistance)) )
  # Numerical vector representing all horizontal tags positioning to be used
  # to correctly locate them on the plot.
  
  text(x = xAxis, # Horizontal positions of features tags. 
       y = -0.01, # Default vertical position of all features tags.
       labels = colnames(barplotData), # Names of features as plot horizontal
       # axis labels.
       xpd = NA, # allows external plot text.
       cex = cexLabels, # Factor that alters labels text size.
       srt = tiltAngle, # Horizontal axis labels tilt angle.
       adj = adjText # Numerical vector that slightly displaces horizontal axis
       # labels.
       ) # plots horizontal axis labels.
  
  legend(x = "topright", # Legend default position on the plot.
         legend = trueClassesOrder, # Classes tags as the legend text.
         col = "black", # Color of legend symbols borders.
         pch = 22, # Type of legend symbol.
         pt.bg = classesColors, # Character vector of ordered colors tags.
         cex = cexLeg, # Factor that determines legend size.
         inset = legInset # Numerical vector that displaces legend across the
         # horizontal and vertical axes.
         ) # plots image legend.
  
  return(list("barplotData" = barplotData, # outputs reordered importance array.
              "percentTable" = percentTableValues) # outputs relative feature
         # percentage importance results.
         )
  
}

## 11) heat:

# This function requires the "gplots" package. An input containing only the
# data in the form of a matrix or data frame is needed, as well as the distance
# calculation method as a string. The data only have to contain the samples and
# feature levels in their columns and rows, respectively. The possible argument
# inputs for distance calculations are: "euclidean", "manhattan" and
# "correlation". The key index shows the deviation from a mean value in terms
# of Z-score (number of standard deviations) calculated using the feature
# levels from all samples. In each axis, a dendrogram can be plotted. While
# the used distance function can be specified, clustering method cannot.
# "heat" function does require NA values omission.

heat <- function(dataF, # Data set having the samples distributed across
                 # columns and its features being its rows.
                 distance = "correlation", # Type of distance measurement
                 # calculation that is going to be run to generate the matrix
                 # dendrograms. Possible distance measurement methods are:
                 # "euclidean", "manhattan" and "correlation".
                 xcex = 0.75, # Factor that alters horizontal axis labels size.
                 ycex = 0.75, # Factor that alters vertical axis labels size.
                 ylab = TRUE, # indicates if vertical axis labels are wanted to
                 # be plotted.
                 rowD = TRUE, # indicates if features clustering dendrogram is
                 # wanted to be calculated.
                 colD = TRUE, # indicates if samples clustering dendrogram is
                 # wanted to be calculated.
                 colSet = "greenred", # indicates the color palette that is
                 # going to be used. From smaller values to greater ones, the
                 # two colors shown in the string represent the limits of a
                 # colors range that is going to be used in the color code
                 # depiction. There are three scripted possibilities:
                 # "greenred", "redgreen" or "bluered".
                 rotateXLabs = NULL, # uses a numerical value to apply a degree
                 # of rotation to the horizontal axis plotted labels.
                 rotateYLabs = NULL, # uses a numerical value to apply a degree
                 # of rotation to the vertical axis plotted labels.
                 ADJCOL = c(NA, 0), # Numerical vector which can be input with
                 # values from -1 to 1, that displaces horizontal axis labels
                 # from left to right and up and down.
                 ADJROW = c(0, NA) # Numerical vector which can be input with
                 # values from -1 to 1, that displaces vertical axis labels
                 # from left to right and up and down.
                 ){
  
  if(!("package:gplots" %in% base::search())){ # If the "gplots" package has
    # not been loaded yet, it loads it into R.
    library(gplots) # attaches "gplots" package.
  }
  
  if(colSet == "greenred"){ # If the color palette has been determined to range
    # from green to red through "colSet" argument, this conditional is executed.
    colsUsed <- colorRampPalette(colors = c("green", "black", "red"))(75)
    # creates a character vector of the scripted length of color names starting
    # from the one indicated first and gradually passing from the current color
    # to the next one in the vector.
  }
  if(colSet == "bluered"){ # If the color palette has been determined to range
    # from blue to red through "colSet" argument, this conditional is executed.
    colsUsed <- colorRampPalette(colors = c("blue", "black", "red"))(75)
    # creates a character vector of the scripted length of color names starting
    # from the one indicated first and gradually passing from the current color
    # to the next one in the vector.
  }
  if(colSet == "redgreen"){ # If the color palette has been determined to range
    # from red to green through "colSet" argument, this conditional is executed.
    colsUsed <- colorRampPalette(colors = c("red", "black", "green"))(75)
    # creates a character vector of the scripted length of color names starting
    # from the one indicated first and gradually passing from the current color
    # to the next one in the vector.
  }
  
  if(ylab){ # If it has been indicated that vertical axis labels are wanted to
    # be plotted, this part of the conditional will be run.
    if(distance == "euclidean" | distance == "manhattan"){ # If the
      # hierarchical dendrogram clustering distance measurement has been
      # indicated to either be based on the computation of euclidean distance
      # or Manhattan distance, this part of the code will be run.
      heatmap.2(x = dataF, # Data set.
                distfun = function(x) {dist(x, method = distance)}, # describes
                # distance calculation function.
                trace = "none", # determines if lines are wanted to be drawn
                # as to define an array.
                scale = "row", # indicates which data set dimension has to be
                # chosen to compute represented Z-scores.
                density.info = "none", # indicates through a string if extra
                # information is wanted on the color key that is attached to
                # the heat map image.
                key.title = "", # String that is incorporated as the color
                # key title.
                cexRow = ycex, # Factor that controls vertical axis labels
                # size.
                cexCol = xcex, # Factor that controls horizontal axis labels
                # size.
                Rowv = rowD, # Boolean value that determines if a dendrogram
                # clustering rows is wanted to be plotted.
                Colv = colD, # Boolean value that determines if a dendrogram
                # clustering columns is wanted to be plotted.
                col = colsUsed, # Set of colors that are going to be used on
                # the heat map to represent the range of values.
                srtCol = rotateXLabs, # Number indicating a degree of rotation
                # that has to be applied on the horizontal axis labels to
                # plot them in the desired orientation.
                srtRow = rotateYLabs, # Number indicating a degree of rotation
                # that has to be applied on the vertical axis labels to plot
                # them in the desired orientation.
                adjCol = ADJCOL, # Numerical vector which can be input with
                # values from -1 to 1, that displaces horizontal axis labels
                # from left to right and up and down. 
                adjRow = ADJROW # Numerical vector which can be input with
                # values from -1 to 1, that displaces vertical axis labels
                # from left to right and up and down.
                )
    }
    if(distance == "correlation"){ # If the hierarchical dendrogram clustering
      # distance measurement has been indicated to be based on the computation
      # of the Pearson's correlation coefficient as distance measurement, this
      # code chunk will be run.
      heatmap.2(x = dataF, # Data set.
                distfun = function(x)
                  as.dist((1 - cor(t(x), use = "pairwise.complete.obs"))/2),
                # describes distance calculation function.
                trace = "none", # determines if lines are wanted to be drawn
                # as to define an array.
                scale = "row", # indicates which data set dimension has to be
                # chosen to compute represented Z-scores.
                density.info = "none", # indicates through a string if extra
                # information is wanted on the color key that is attached to
                # the heat map image.
                key.title = "", # String that is incorporated as the color
                # key title.
                cexRow = ycex, # Factor that controls vertical axis labels
                # size.
                cexCol = xcex, # Factor that controls horizontal axis labels
                # size.
                Rowv = rowD, # Boolean value that determines if a dendrogram
                # clustering rows is wanted to be plotted.
                Colv = colD, # Boolean value that determines if a dendrogram
                # clustering columns is wanted to be plotted.
                col = colsUsed, # Set of colors that are going to be used on
                # the heat map to represent the range of values.
                srtCol = rotateXLabs, # Number indicating a degree of rotation
                # that has to be applied on the horizontal axis labels to
                # plot them in the desired orientation.
                srtRow = rotateYLabs, # Number indicating a degree of rotation
                # that has to be applied on the vertical axis labels to plot
                # them in the desired orientation.
                adjCol = ADJCOL, # Numerical vector which can be input with
                # values from -1 to 1, that displaces horizontal axis labels
                # from left to right and up and down. 
                adjRow = ADJROW # Numerical vector which can be input with
                # values from -1 to 1, that displaces vertical axis labels
                # from left to right and up and down.
                )
    }
  }else{ # If it has been indicated that vertical axis labels are wanted not to
    # be plotted, this part of the conditional will be run.
    if(distance == "euclidean" | distance == "manhattan"){ # If the
      # hierarchical dendrogram clustering distance measurement has been
      # indicated to either be based on the computation of euclidean distance
      # or manhattan distance, this part of the code will be run.
      
      # With respect to previous code chunks, these next ones only differ in
      # the presence of the argument "labRow". As it is indicated in R help
      # page, labRow represents: "character vector with row labels to use".
      # For understanding the meaning of the rest of the code, check the
      # code lines above.
      
      heatmap.2(x = dataF,
                distfun = function(x) {dist(x, method = distance)},
                trace = "none",
                scale = "row",
                density.info = "none",
                key.title = "",
                cexRow = ycex,
                cexCol = xcex,
                labRow = ylab, # Character vector with row labels to use.
                Rowv = rowD,
                Colv = colD,
                col = colsUsed,
                srtCol = rotateXLabs,
                srtRow = rotateYLabs,
                adjCol = ADJCOL,
                adjRow = ADJROW)
    }
    if(distance == "correlation"){ # If the hierarchical dendrogram clustering
      # distance measurement has been indicated to be based on the computation
      # of the Pearson's correlation coefficient as distance measurement, this
      # code chunk will be run.
      heatmap.2(x = dataF,
                distfun = function(x)
                  as.dist((1 - cor(t(x), use = "pairwise.complete.obs"))/2),
                trace = "none",
                scale = "row",
                density.info = "none",
                key.title = "",
                cexRow = ycex,
                cexCol = xcex,
                labRow = ylab, # Character vector with row labels to use.
                Rowv = rowD,
                Colv = colD,
                col = colsUsed,
                srtCol = rotateXLabs,
                srtRow = rotateYLabs,
                adjCol = ADJCOL,
                adjRow = ADJROW)
    }
  }
}

## 12) plot3D:

# This function was created to generate 3D scatter plots in R. A data set
# having three columns is required, so the data from each column represents
# the data related to one dimension.

plot3D <- function(data3D, # Data arranged in a three-column matrix or data
                   # frame.
                   classVector, # A character vector having the samples classes
                   # tags ordered in such a way so they correspond to the
                   # data array or data frame ordered rows. 
                   classColors, # A character vector with as many elements as
                   # unique classes are. Each color will be orderly assigned to
                   # a class.
                   main3D = "3D plot", # Plot title.
                   cex3D = 1.2, # Factor that controls point size on the plot.
                   cexLegend = 1.5, # Factor that controls plot legend size.
                   xAxisLimits = NULL, # Numerical vector to indicate the
                   # lowest and the highest value associated the "x" axis
                   # dimension.
                   yAxisLimits = NULL, # Numerical vector to indicate the
                   # lowest and the highest value associated the "y" axis
                   # dimension.
                   zAxisLimits = NULL, # Numerical vector to indicate the
                   # lowest and the highest value associated the "z" axis
                   # dimension.
                   par3D = 1 # Factor that controls numerical and character
                   # labels size from the axes.
                   ){
  
  if(!("package:rgl" %in% base::search())){ # Conditional that is run if the
    # "rgl" package is not loaded.
    library(rgl) # loads "rgl" package.
  }
  
  classesOnlyOnce <- unique(classVector) # stores the classes tags without
  # repetition.
  namedColors <- classColors # creates a replica of the classes colors
  # character vector.
  names(namedColors) <- classesOnlyOnce # names colors vector using classes
  # vector with no repeated elements.
  completeColorVector <- namedColors[classVector] # creates a color vector
  # which length equals the number of rows the input data set has. This vector
  # corresponds to the color that should represent each sample considering
  # the class it belongs to.
  
  if(is.null(xAxisLimits)){ # If the "x" axis limits have not been indicated
    # when executing the function, this code chunk will be executed.
    xDif <- pretty(data3D[,1])[2] - pretty(data3D[,1])[1] # "pretty" function
    # creates a numerical vector if two values are conceded as input. The
    # result is a series of consecutive values that are equally distanced
    # between each other and constitute a range in which the original data can
    # be found. This constant difference between consecutive values is
    # calculated by subtracting the second and first vector elements.
    xRangeMin <- min(data3D[,1]) - xDif # subtracts previous difference to
    # "x" axis data minimal value.
    xRangeMax <- max(data3D[,1]) + xDif # adds previous difference to "x" axis
    # data maximal value.
    xRange <- range(pretty(c(xRangeMin, xRangeMax))) # creates a new range of
    # values from inputting two last computed number into the "pretty"
    # function. After that, the minimum and maximum values are stored as "x"
    # axis limits.
  }else{ # If the "x" axis limits have already been indicated in the
    # "xAxisLimits" function argument, this part of the conditional will be
    # run:
    xRange <- c(xAxisLimits[1], xAxisLimits[2]) # uses input values as "x" axis
    # limits.
  }
  
  # The same code structure is applied to the "y" axis definition and the "z"
  # axis definition.
  
  if(is.null(yAxisLimits)){
    yDif <- pretty(data3D[,2])[2] - pretty(data3D[,2])[1]
    yRangeMin <- min(data3D[,2]) - yDif
    yRangeMax <- max(data3D[,2]) + yDif
    yRange <- range(pretty(c(yRangeMin, yRangeMax)))
  }else{
    yRange <- c(yAxisLimits[1], yAxisLimits[2])
  }
  
  if(is.null(zAxisLimits)){
    zDif <- pretty(data3D[,3])[2] - pretty(data3D[,3])[1]
    zRangeMin <- min(data3D[,3]) - zDif
    zRangeMax <- max(data3D[,3]) + zDif
    zRange <- range(pretty(c(zRangeMin, zRangeMax)))
  }else{
    zRange <- c(zAxisLimits[1], zAxisLimits[2])
  }
  
  if(is.null(dimensionNames[1])){
    dimensionNames[1] <- "Dimension 1"
  }
  if(is.null(dimensionNames[2])){
    dimensionNames[2] <- "Dimension 2"
  }
  if(is.null(dimensionNames[3])){
    dimensionNames[3] <- "Dimension 3"
  }
  
  par3d(cex = par3D) # modifies 3D environment to make the text from axes
  # labels adapted to a wanted size based on a given number.
  
  plot3d(x = data3D[,1], # Data to be plotted on "x" axis.
         y = data3D[,2], # Data to be plotted on "y" axis.
         z = data3D[,3], # Data to be plotted on "z" axis.
         main = main3D, # Plot title.
         xlab = dimensionNames[1], # "x" axis name label.
         ylab = dimensionNames[2], # "y" axis name label.
         zlab = dimensionNames[3], # "z" axis name label.
         type = "s", # Type of point to be represented. "s" stands for spheres.
         size = cex3D, # Factor that modulates point size.
         col = completeColorVector, # Complete character vector indicating
         # the color of each sample in order.
         xlim = xRange, # "x" axis limits.
         ylim = yRange, # "y" axis limits.
         zlim = zRange # "z" axis limits.
         ) # creates 3D plot.
  
  legend3d(x = "topright",
           legend = classesOnlyOnce,
           col = classColors,
           pch = 19,
           cex = cexLegend,
           inset = c(0.05))
  
}

## 13) PCAGif:

# It creates a GIF file from a 3D interactive plot. Using this function, it can
# be indicating in which sense the plot is going to rotate to create a movie of
# it.

PCAGif <- function(sense = c(0, 0, 1), # A numerical vector with three
                   # elements, each one of them representing an axis. A zero
                   # represents no move generated from rotating that axis,
                   # while a -1 or 1 represents movement, in one sense or the
                   # other. Each position in the vector refers to one of the
                   # dimensions.
                   nameMe = "PCAGif" # names GIF file (GIF extension is given
                   # to the file name by default).
                   ){
  
  if(!("package:rgl" %in% base::search())){ # If the "rgl" library has not been
    # loaded, this conditional is executed.
    library(rgl) # loads "rgl" library.
  }
  
  movie3d(spin3d(axis = sense, # indicates rotation sense.
                 rpm = 4 # Number of revolutions per minute.
                 ), # makes plot rotate.
          duration = 15, # Number of seconds the movie will last (4 rpm is
          # equivalent to 1 revolution in 1/4 minutes, which is 15 seconds.
          # This means that the movie is programmed to capture a single plot
          # spin).
          dir = "./", # Saving path.
          movie = nameMe # String to name the GIF file without explicitly
          # mentioning its format.
          )

}

