# This script handles a data table of mass spectrometry results to eventually
# find biomarkers needed to differentiate samples coming from different
# patients showing different diseases but with very similar symptoms. Through
# a proteomics analysis, differences can be found and, at the same time,
# patients from the same class due to having the same illness can be grouped
# together and characterized.

# This script has been plainly divided into 6 sequential sections. In each one,
# a part of the analysis is executed. These 6 parts have been named as follows:

# 1. Determining working directory
# 2. Loading necessary scripted functions and libraries in current computer
# 3. Loading file information
# 4. Data standardization
# 5. Application of machine learning algorithms
# 6. Extra images output

################################################################################

# 1) Determining working directory

workDirPath <- "Z:/USUARIOS/JVINDEL/laura/cholestasis20230711/"
# establishes complete working path name.
setwd(workDirPath) # sets it.
dir.create(path = "rda") # creates another folder into current working
# directory.

# The computer working path. It can be changed depending on the system where
# the script is being executed.

# 2) Loading necessary scripted functions and libraries in current computer

cholestasisProteomeFunctionsFilePath <-
  "C:/Users/Centaurus2/Documents/cholestasisComplete_functions.R"
# All scripted functions to perform this analysis can be found in the
# previous file. The directory can be changed if needed.
source(cholestasisProteomeFunctionsFilePath)
# It runs the script that contains the scripted functions for this analysis.

allLibraries <- c("limma", "proBatch",
                  "RColorBrewer", "gplots", "DescTools", "xlsx", "randomForest",
                  "plot.matrix", "caret", "MASS", "e1071", "xgboost",
                  "colorspace")
# Name of the libraries that are going to be required for this analysis. To
# install or attach them, the "checkPackage" function is going to be used. It'll
# install the libraries, or just attach them if they are already available on
# the current computer.
allRepos <- c(rep("Bioconductor", 2),
              rep("CRAN", 11)
              ) # Character vector that stores the repository each library
# comes from. This is necessary since the "checkPackage" function works using
# the name of the library and also the name of the repository it was
# deposited in. Even though the repository name is only strictly necessary
# when it has to be installed, it allows this function to install and attach
# all the libraries if needed just by only running it in a loop.

for(i in 1:length(allLibraries)){ # It iteratively installs or attaches each
  # library.
  checkPackage(packageLabel = allLibraries[i], # selected library.
               repository = allRepos[i] # repository the library belongs to.
  )
}

# 3) Loading file information

# Handling full data table:

proteomeDataPath <- paste("originalFiles/",
                          "20220307_TMT_cholestasis_all_high_FDR_proteins.csv",
                          sep = "")
# Complete proteome data file path in current computer.

proteomeData <- read.csv2(proteomeDataPath) # loads data file.
nrow(proteomeData) # Original set of proteins: 7161 proteins.
proteomeData <- proteomeData[proteomeData[,3] == "Master Protein",]
# Only proteins detected as "master proteins" will be considered in the
# analysis.
nrow(proteomeData) # Number of "master proteins": 6648.

proteomeData[grep(pattern = "-Cont)",
                  x = proteomeData[, "Description"],
                  fixed = TRUE),]
proteomeData <- proteomeData[-grep(pattern = "-Cont)",
                                   x = proteomeData[, "Description"],
                                   fixed = TRUE),]
# Protein contaminants will be deleted from the data.
nrow(proteomeData)
# number of "master" proteins that are not contaminants: 6641.
proteomeData[grep(pattern = "MW-Marker",
                  x = proteomeData[, "Description"],
                  fixed = TRUE),]
proteomeData <- proteomeData[-grep(pattern = "MW-Marker",
                                   x = proteomeData[, "Description"],
                                   fixed = TRUE),]
# Also, molecular biomarkers should not be considered into the analysis and
# therefore can be deleted from the data table.
nrow(proteomeData) # number of "master" proteins that are not either
# contaminants or molecular biomarkers: 6640.

# Handling description data table:

proteomeDescription <- proteomeData[,1:2] # retrieves only accession codes and
# proteins description as a data frame.
proteomeDescription$gene.ID <- geneIDextraction(proteomeDescription[,2])
# retrieves gene codes from proteins description and merges protein accession
# codes, protein descriptions and gene code names into one single data frame.
colnames(proteomeDescription) <- c("accession", "description", "gene.ID")
# customs previous data frame columns names.
rownames(proteomeDescription) <- proteomeDescription[,1] # names rows after
# UniProt accession codes. This will make row indexation easier.

save(x = proteomeDescription, file = "rda/proteomeDescription.rda")
# saves proteins description table as a RDA object in a specific path.

# Abundances data table:

proteomeAbundances <- proteomeData[, 4:ncol(proteomeData)] # separates samples
# data columns from the rest of the data.
rownames(proteomeAbundances) <- rownames(proteomeDescription)
# names protein abundances data array rows after protein description table
# rows.

# To change the name of samples that are shown as column names in the data
# array, a copy of them as a character vector is going to be created and
# some substrings in this vector will be replaced, to make nomenclature
# more abbreviated and comprehensible.

colsCopy <- colnames(proteomeAbundances) # creates a copy of the character
# vector of samples names.

patterns <- c("Abundance..", # 1
              "..Sample..", # 2
              "Alagille.Syndrome", # 3
              "Biliary.Atresia", # 4
              "..Control..IS", # 5
              "..", # 6
              "DAAT" # 7
              ) # patterns that are wanted to be replaced in order.
replacements <- c("", # 1
                  ".", # 2
                  "Alagille", # 3
                  "Atresia", # 4
                  ".IS", # 5
                  ".", # 6
                  "AATD" # 7
                  ) # pattern replacements in order.

for(i in 1:length(patterns)){ # loop that iterates over all patterns and
  # their corresponding replacements.
  
  colsCopy <- gsub(pattern = patterns[i],
                   replacement = replacements[i],
                   x = colsCopy,
                   fixed = TRUE) # replaces each pattern by its replacement in
                                 # order.
  
}

colnames(proteomeAbundances) <- colsCopy # substitutes default samples names
# by their shortened name version.

proteomeAbundances <- as.matrix(proteomeAbundances) # protein abundance data
# is converted from being a data frame into an array.
nrow(proteomeAbundances) # checks number of rows in current array: 6640.
proteomeAbundances <- na.omit(proteomeAbundances) # deletes those rows in
# the data array that have at least one missing value.

# To process the data, it has been decided that those proteins that present
# at least one missing value will be erased from the data set. Some approaches
# use data imputation, but due to the great number of features and the small
# number of samples, just omitting rows that have some missing value seems
# a sensible approach. Other observation that is worth noticing is that missing
# values tend to be grouped by samples batch.

nrow(proteomeAbundances) # checks number of rows in current array: 5487.
# At least 6640 - 5487 = 1153 proteins had a single missing value.

sum(proteomeAbundances < 1) # calculates the number of data values that are
# less than 1 as raw protein abundance measurement. There are 5 data that are
# less than 1.
sum(proteomeAbundances < 1) == sum(proteomeAbundances == 0) # checks if those
# data specifically equal zero, which is true.

proteomeAbundances[proteomeAbundances == 0] <- NA # Those data which
# equal zero are considered missing values, so they are replaced by "Not
# Available" values.
proteomeAbundances <- na.omit(proteomeAbundances) # New rows with missing values
# are removed.
nrow(proteomeAbundances) # current number of rows: 5484.

# Now, data array columns are going to be reorganized by groups, in a way that
# samples are presented grouped by classes and classes are arranged in the
# following order: control group, Alagille syndrome, biliary atresia,
# alpha-1-antitrypsin deficit (AATD) and progressive familiar intrahepatic
# cholestasis (PFIC). To do so, multiple functions to get indexes by substring
# presence and setting differences between other R objects from the same class
# are going to be used in the next code lines:

proteomeAbundances <- proteomeAbundances[,c(grep(pattern = "Control",
                                                 x = colnames(proteomeAbundances)),
                                            setdiff(x = setdiff(x = 1:ncol(proteomeAbundances),
                                                                y = grep(pattern = "Control",
                                                                         x = colnames(proteomeAbundances))),
                                                    y = grep(pattern = "IS",
                                                             x = colnames(proteomeAbundances))),
                                            grep(pattern = "IS",
                                                 x = colnames(proteomeAbundances))
                                            )]

# Basically, considering samples classes known distribution across data table
# columns, control group samples indexes are obtained to enumerate them as the
# first data columns. After that, indexes from all classes except the control
# group and the internal standards group are retrieved. Since these samples are
# ordered as desired with respect to previous explained samples ordering,
# samples groups are added like that to the new data table and finally,
# internal standards samples are included as the last columns of the data table.

save(x = proteomeAbundances, file = "rda/proteomeAbundances.rda") # saves
# raw ordered data result. 

# 4. Data standardization

# Logarithmic transformation is common when operating with transcriptomics,
# proteomics or metabolomics data. Logarithmic transformation contributes to
# outlier smoothing and data adopting a normal distribution, which is useful
# when performing parametric statistical tests.

log2Proteome <- log2(proteomeAbundances) # Data logarithmic transformation is
# calculated.

# 4 different batches constitute the data, each one having its own internal
# standard. Using it, batch effect will be removed to be able to compare all
# data between batches.

# Batch belonging is represented by a tag in the samples names. The samples
# naming is due to data having being analysed by Proteome Discoverer software.
# Information from each batch is located in a file ("F"). For all batches in an
# experiment, there will be a number starting from 1 that will also identify
# them. Using regular expressions, batch tag can be obtained from samples names:

batches <- regexpr(pattern = "F[0-9]{1}", # regular expression
                   text = colnames(proteomeAbundances) # text onto which
                   # patterns will be searched.
                   ) # gets patterns character positioning indices.
batches <- substr(x = colnames(proteomeAbundances),
                  start = batches,
                  stop = batches + 
                    attr(x = batches,
                         which = names(attributes(batches))[1]) - 1)
# gets substrings from dividing strings into chunks by indicating character
# positions.

classesVector <- unlist(lapply(X = strsplit(x = colnames(proteomeAbundances),
                                            split = ".",
                                            fixed = TRUE),
                               FUN = function(x) x[length(x)]))
# retrieves what will be the class tags from samples names.

uniqueClasses <- unique(classesVector) # deletes all repeated elements in a 
# vector so all elements only appear once.

batchRemoved <- log2Proteome # creates a copy of the logarithmically
# transformed data to use it as a scaffold to remove batch effect.

ISBool <- grep(pattern = "IS",
               x = colnames(batchRemoved),
               fixed = TRUE) # generates logical vector to only index
# data columns that are marked as internal standards.

for(i in 1:length(unique(batches))){ # loop with as many iterations as
  # batches are.
  
  ibatch <- unique(batches)[i] # retrieves specific batch tag.
  ibatchBool <- grep(pattern = ibatch,
                     x = colnames(batchRemoved),
                     fixed = TRUE) # generates logical vector to only index
  # data columns that belong to the specific batch.
  
  iISBool <- intersect(x = ibatchBool, y = ISBool) # intersects samples indices
  # from a specific batch and those indices of internal standards to obtain
  # only the sample index that represents the internal standard of a batch. 
  iIScolName <- colnames(batchRemoved)[iISBool] # the selected batch internal
  # standard sample name is recorded.
  iISdata <- batchRemoved[,iISBool] # saves the data column of the specific
  # batch internal standard.
  ibatchData <- batchRemoved[,ibatchBool] # separates data samples from a 
  # specific batch from the rest. This constitutes an array made of a subset of
  # samples.
  ibatchData <- ibatchData[, -grep(pattern = "IS",
                                   x = colnames(ibatchData),
                                   fixed = TRUE)] # deletes the internal
  # standard from all the separated samples data of a batch.
  
  for(j in 1:ncol(ibatchData)){ # loop that iterates over all previous batch
    # samples.
    
    jcol <- ibatchData[, j] - iISdata # subtracts the batch internal standard
    # from all samples of the same batch, one by one. This removes the batch
    # effect from the samples.
    jcolName <- colnames(ibatchData)[j] # retrieves sample name. By doing so,
    # batch removal results can be integrated back into the original data array.
    batchRemoved[, jcolName] <- jcol # uses sample name or column name to
    # introduce the sample data without batch effect back into the array
    # with all the samples.
    
  }
  
  iISdata <- iISdata - iISdata # controls that internal standard data have been
  # used to correctly remove the batch effect from samples and finally
  # subtracts itself to equal zero for all of its values. By doing so on all
  # internal standard data, the result has to show that all internal standards
  # are coincident by their values.
  batchRemoved[, iIScolName] <- iISdata # returns updated internal standard data
  # into the original data matrix.
  
}

save(x = batchRemoved, file = "rda/batchRemoved.rda") # saves data after
# removing all batch effect from them.

firstColorSet <- c(RColorBrewer::brewer.pal(n = 8, name = "Set1")[c(3,2,5,1)],
                   "gold", "darkgrey") # selects a set of colors to represent
# each of the sample classes available in this data set.

samplesColorVector <- c() # Empty vector that will be used to assign one color
# to each sample. For all samples from a class, different color tone similar
# between them will be used.
colorChangeFactor <- 0.15 # A value used to lighten or darken colors associated
# to samples from the same class, in order to be able to group them by color
# similarity and, at the same time, to be able to visually differentiate every
# sample from the rest.

for(i in 1:length(firstColorSet)){ # This loop will be used to generate all
  # colors associated to a class by assigning one color to each sample. In
  # each round, all colors related to a class will be chosen.
  
  icolor <- firstColorSet[i] # A base color for the class is selected.
  iclass <- unique(classesVector)[i] # A class is selected.
  ncolors <- table(classesVector)[iclass] # The number of samples belonging to
  # the class is obtained.
  ireps <- floor(ncolors/2) # The number of samples of this class is divided
  # by 2 and the result is rounded to the smallest number.
  ilightColors <- lighten(col = icolor, amount = rep(x = colorChangeFactor, ireps) * 1:ireps)
  # Using the base class color and the number of samples per class halved,
  # an equal number of lightened colors is obtained. This set of colors will be
  # assigned to a subset of the samples from a specific class.
  ilightColors <- rev(ilightColors) # The character vector containing all
  # lightened colors is reversed considering colors order.
  if(ireps %% 2 == 0){ # Depending on an even or odd number of samples in
    # classes, next a set of darkened colors will be produced in the same way
    # lightened color were generated.
    ireps <- ireps + 1
    idarkColors <- darken(col = icolor, amount = rep(x = colorChangeFactor, ireps) * 1:ireps)
  }else{
    idarkColors <- darken(col = icolor, amount = rep(x = colorChangeFactor, ireps) * 1:ireps)
  }
  samplesColorVector <- append(x = samplesColorVector,
                               values = c(ilightColors, icolor, idarkColors))
  # All colors are added to the original empty vector created to hold a color
  # associated to each sample.
  
}

medianNorm <- proBatch::normalize_data_dm(data_matrix = batchRemoved,
                                          normalize_func = c("medianCentering"))
# Data matrix is standardized using "proBatch" package normalization function.
# Chosen data normalization is median normalization.
medianNorm <- medianNorm[, -grep(pattern = "IS", x = colnames(medianNorm))]
# deletes internal standard data from normalized data matrix.

# Here below, a set of different figures representing samples data at the
# different transformation and standardization stages will be plotted. Used
# data arrays will be: original data matrix, logarithmically transformed data
# set, removed batch effect data and finally median-standardized data array.
# For each array, PCA, data distribution curves and boxplots will be plotted
# in R. This will help to visualize data structures across data processing
# stages.

for(i in 1:length(list(proteomeAbundances, log2Proteome, batchRemoved, medianNorm))){
  
  # This loop iterates over multiple data matrices after having altered the data
  # in different steps in order to process the data. By transforming and
  # standardizing data, it is expected to derive from them the real knowledge
  # that is condensed in them, ending with results that would not represent
  # false positives or false negatives.
  
  idata <- list(proteomeAbundances, log2Proteome, batchRemoved, medianNorm)[[i]]
  # A data matrix is selected from the processed data compilation.
  iDensityCurvesTitle <- c("Samples histogram density curves (no log. transf.; no batch correction; no norm.)",
                           "Samples histogram density curves (log. transf.; no batch correction; no norm.)",
                           "Samples histogram density curves (log. transf.; IS correction; no norm.)",
                           "Samples histogram density curves (log. transf.; IS correction; median norm.)")[i]
  # A data density curve title is chosen from the manually written list.
  iPCAtitle <- c("PCA (proteome data; all proteins) (no log. transf.; no batch correction; no norm.)",
                 "PCA (proteome data; all proteins) (log. transf.; no batch correction; no norm.)",
                 "PCA (proteome data; all proteins) (log. transf.; IS correction; no norm.)",
                 "PCA (proteome data; all proteins) (log. transf.; IS correction; median norm.)")[i]
  # A PCA title is chosen from the manually written list.
  idensityall <- density(idata) # Data density curve is generated from a
  # data matrix.
  ixlimits <- range(pretty(idensityall$x)) # ranges data spreading to make
  # an aesthetically pleasing plot by limiting horizontal axis.
  iylimits <- c(0, 0) # Vertical axis range variable shell. It'll be updated
  # with the best values to plot an aesthetic image.
  
  for(j in 1:ncol(idata)){ # This loop iterates over all data matrix samples
    # to find the one with the greatest associated density curve results.
    # By considering these higher values, the density curve plot vertical axis
    # will be correctly adjusted.
    
    jylimits <- range(pretty(density(idata[,j])$y)) # ranges sample density
    # curve values to establish a maximum that allows correct plotting.
    if(jylimits[2] > iylimits[2]){ # Vertical axis range iteratively updates
      # if a sample shows a higher density curve maximum.
      iylimits <- jylimits
    }
    
  }
  
  for(j in 1:ncol(idata)){ # Loop that iterates over all samples to plot
    # their density curves.
    
    jdensity <- density(x = idata[,j]) # returns density curve data.
    if(j == 1){
      plot(jdensity,
           lty = 1,
           lwd = 2,
           col = samplesColorVector[j],
           main = iDensityCurvesTitle,
           xlab = "Data",
           ylab = "Relative frequency",
           xlim = ixlimits,
           ylim = iylimits
           ) # plots data density results as a curve. Note that each sample
      # will have a corresponding assigned color considering its class and
      # a slightly lightened or darkened tone with respect to its class
      # original color.
    }else{
      points(jdensity,
             type = "l",
             lty = 1,
             lwd = 2,
             col = samplesColorVector[j]) # plots the rest of the density
      # curves over the first curve plot.
    }
    
  }
  
  legend(x = "topright",
         legend = colnames(idata),
         col = samplesColorVector,
         lty = 1,
         lwd = 2,
         cex = 0.565) # adds a legend to density curves plot.
  
  iboxTitle <- c("Raw values boxplots of samples",
                 "log2 transformed values boxplots of samples",
                 "No batch effect values boxplots of samples",
                 "Complete preprocessed values boxplots of samples")[i]
  # orderly chooses a boxplot title associated to its corresponding data array.
  iyBoxLab <- c("Raw values",
                "log2 transformed values",
                "No batch effect values",
                "Complete preprocessed values")[i]
  # Each boxplot has its adequate vertical axis label depending on the data
  # transformation or standardization stage.
  
  boxplot(idata,
          col = samplesColorVector,
          pch = 21,
          bg = samplesColorVector,
          main = iboxTitle,
          ylab = iyBoxLab,
          xaxt = "n") # generates boxplot image to represent data values per
  # sample.
  
  iy <- c(-17000, -0.2, -5.6, -5.6)[i] # scripted text vertical position to be
  # used to include text as sample name labels.
  
  axis(side = 1, at = 1:ncol(idata), labels = rep("", ncol(idata)))
  # plots horizontal axis label marks for all boxplots.
  text(x = 1:ncol(idata),
       y = iy,
       labels = colnames(idata),
       cex = 0.8,
       srt = 45,
       adj = 1,
       xpd = TRUE) # plots samples names, one for each corresponding boxplot.
  
  if(i == 4){ # If during the current iteration the normalized data set by
    # median normalization is being used, the last group in order, which is the
    # internal standards group, will be deleted and will not be considered in
    # the PCA plot.
    iclassesForPCA <- uniqueClasses[-length(uniqueClasses)]
    icolorsForPCA <- firstColorSet[-length(firstColorSet)] # It also considers
    # internal standard group color deleting from the group colors vector.
  }else{
    iclassesForPCA <- uniqueClasses
    icolorsForPCA <- firstColorSet
  }
  ilegPosPCA <- c("bottomright",
                  "bottomright",
                  "topleft",
                  "bottomleft")[i]
  # scriptedly chooses legend position for PCA plot.
  
  PCAPlot(PCAObject = PCACal(scaledData = idata, wanna.scale = FALSE),
          is.chosen = TRUE,
          firstPC = 1,
          secPC = 2,
          groups = iclassesForPCA,
          colVec = icolorsForPCA,
          legPos = ilegPosPCA,
          main2 = iPCAtitle,
          cexLeg = 1.3,
          cexPoints = 1.5) # calls a user-defined function to make a PCA plot.
  # It considers the different classes in the data set to adequately color
  # samples, adds a title and a legend, and plotted PC can be selected.
  
}

save(x = medianNorm, file = "rda/medianNorm.rda") # saves data after applying
# median normalization as an RDA file.
save(x = classesVector, file = "rda/classesVector.rda") # saves ordered class
# belonging vector as an RDA file.

noISclassesVector <- unlist(lapply(X = strsplit(x = colnames(medianNorm),
                                                split = ".",
                                                fixed = TRUE),
                                   FUN = function(x) x[length(x)]))
# extracts class belonging from median-normalized data array column names,
# and so it does not consider internal standard samples.

save(x = noISclassesVector, file = "rda/noISclassesVector.rda") # saves
# previous classes vector as an RDA file in current working directory.

# 5) Application of machine learning algorithms

load(file = "rda/medianNorm.rda") # loads median-normalized data.
load(file = "rda/noISclassesVector.rda") # loads classes vector data.

cholestasisData <- cbind.data.frame(t(medianNorm),
                                    "patient" = as.factor(noISclassesVector))
# creates a data frame containing full protein data and classes belonging
# data using previous classes vector.
partitionProportion <- 0.7 # Proportion of data to be separated as training set.
set.seed(1) # establishes arbitrary programming seed environment.
partitionIndices <- as.numeric(createDataPartition(y = as.factor(noISclassesVector),
                                                   p = partitionProportion,
                                                   list = FALSE))
# creates a random selection of indices which constitute a subset corresponding
# to the 70% of all data to separate training set from test set,
# which will be used in the next machine learning computational analyses.
cholestasisTrainingData <- cholestasisData[partitionIndices,]
# creates the training set from the full data frame after having randomly
# chosen a subset of samples indices.
cholestasisTestData <- cholestasisData[-partitionIndices,]
# creates test set from remaining sample indices.
cholestasisTrainingClasses <- noISclassesVector[partitionIndices]
# separates training set class tags from full classes vector.
cholestasisTestClasses <- noISclassesVector[-partitionIndices]
# separates test set class tags from full classes vector.

# First, building a classifier using only the training set to select those
# features that seem important for classification will be attempted. Based
# on the results that follow, considering alternative classifier building
# methods would be necessary to extract as much knowledge as possible from
# available data.

set.seed(1) # establishes arbitrary programming seed environment just before
# running the user-defined random forest-based recursive feature elimination
# algorithm.
# The next function has multiple arguments:
fortunaCholestasis <- fortuna(dataset = t(medianNorm[, partitionIndices]),
                              # Samples set from which to analyse feature
                              # relevance.
                              groups = as.factor(cholestasisTrainingClasses),
                              # Class vector as a factor object.
                              treeNumber = 1000, # Number of trees each forest
                              # will have.
                              retentionByProp = FALSE, # indicates if the
                              # number of feature elimination rounds is given
                              # by a factor which is a proportion of the
                              # preserved data per round ("TRUE"), or in
                              # contrast if the number of elimination rounds is
                              # going to be explicitly conceded ("FALSE").
                              retentionProportion = 0.8, # Proportion of the
                              # preserved data per round if "retentionByProp"
                              # argument equals "TRUE".
                              retentionLoops = 30, # Number of elimination
                              # rounds before the entire data set is again
                              # rebuilt to repeat the entire process.
                              chosenIterationsNumber = 30, # Number of times the
                              # entire data set is regenerated after the series
                              # of feature elimination rounds.
                              isBinaryApproach = FALSE, # can have either "TRUE"
                              # or "FALSE" as its value. "isBinaryApproach" was
                              # designed to find in multiclass classification
                              # problems those features that characterized each
                              # class by facing all samples from a class to the
                              # rest of samples that will be included in another
                              # single class. By using
                              # "isBinaryApproach = TRUE", the entire algorithm
                              # will be executed as many times as classes are,
                              # and an importance table result and a classifier
                              # will be generated associated to the
                              # singularization of each class. If
                              # "isBinaryApproach = FALSE", all classes will
                              # be compared at the same time in the
                              # classification algorithm.
                              chosenFeatsClassMethod = "k", # This argument
                              # corresponds to the manner of selecting a subset
                              # of features with the intention of building an
                              # improved classifier based on features
                              # classification performance. This argument can
                              # have as values: "k" (which refers to a given
                              # number of the most important features), "elbow"
                              # (which refers to using the elbow method to
                              # heuristically isolate those features that seem
                              # to be more important than the rest considering
                              # their importance values) or "p" (which refers
                              # to choosing those significant features after
                              # the binomial test significance results).
                              chosenPValue = 0.05, # If
                              # chosenFeatsClassMethod = 'p', the value from
                              # this argument will be considered the
                              # significance threshold.
                              pAdjMethod = "fdr", # Method of p-value
                              # adjustment. Available adjustment methods
                              # correspond to those given by base R "p.adjust"
                              # function.
                              kFeats = 20, # If "chosenFeatsClassMethod = 'k'",
                              # it refers to the number of most important
                              # features to be considered.
                              isElbowPlot = TRUE # If "isElbowPlot = TRUE", a
                              # plot showing importance value decreasing across
                              # features will be generated. It can show that
                              # only a proportion of features is usually useful
                              # when building a classifier, mainly when the
                              # number of features is high.
                              )
# executes random forest feature elimination algorithm considering feature
# importance measurement. This algorithm was primarily based on Diaz-Uriarte's
# discoveries (2006). Later, I realised Boruta algorithm (Kursa, 2010) existed
# and was able to amazingly discriminate between important and useless features
# in classification tasks, relying on chained statistical validations. The
# random forest algorithm presented here has steps that differ from Boruta
# algorithm, such as the reinitialization of the entire data set to repeat all
# processing to make sure that the first deleted features were not erased by
# chance, even though Boruta is faster and provides robust and standardized
# importance results.

save(x = fortunaCholestasis, file = "rda/fortunaCholestasis.rda") # saves
# feature importance results.

dir.create(path = "fortuna")

write.xlsx(x = cbind.data.frame(proteomeDescription[rownames(fortunaCholestasis$importance.matrix$impArr),],
                                fortunaCholestasis$importance.matrix$impArr),
           file = "fortuna/fortunaCholestasis.xlsx",
           sheetName = "Sheet1",
           col.names = TRUE,
           row.names = FALSE,
           append = FALSE)
# uses a built-in function from a package to create XLSX files. The file that
# is going to be created in the form of a table contains protein names,
# descriptions and gene names together with the feature importance results
# from random forest.

# It has to be considered that due to the reduced number of samples,
# conclusions derived from a small set of samples after using machine learning
# methods could be affected by a high variance from the data. Taking into
# account so, the latest importance measurement procedure is going to be
# applied on multiple subsets after different data partitions considering
# sample subsetting. By doing so, average conclusions can be deduced, and these
# results can be compared to those based on individual partitions. This way,
# we can assess if a unique training set is sufficient or maybe by contrast we
# have to consider all samples to obtain more robust answers.

load(file = "rda/fortunaCholestasis.rda") # loads feature importance results
# from the recursive random forest algorithm.

nFirstMostImportant <- 20 # It has been decided that searching among the 20 most
# important features in this classification problem probably harbors a subset
# of features that can provide an appropriate combination to correctly classify
# the data set samples. Considering the small size of the data set, a high
# number of features utilized to supply a classification would be risky to use
# because of overfitting tendencies. Plus, among thousands of proteins from the
# analysis, only a relatively reduced set of them contribute to classification
# with a clearly high importance measurement. The difference of importance
# between features shows that there are only a few among thousands of them
# that are well adapted to this classification problem.
FMIproteins <- rownames(fortunaCholestasis[["importance.matrix"]][["impArr"]])[1:nFirstMostImportant]
# selects the most important proteins after the importance analysis,
# considering the number of most important features conceded in the code line
# immediately above.

load(file = "rda/proteomeDescription.rda") # loads proteins description table.
load(file = "rda/noISclassesVector.rda") # loads classes vector.

seeds <- 100 # Number of partitions that are going to be executed on the data.

allKBestProts <- list() # Empty list that will hold the most important protein
# sets from each data partition.
allPartitionedSamples <- list() # Empty list that will contain the compilation
# of training set samples as indices used after each partition.
dataFramesList <- list() # Empty list that will hold protein importance results
# for each partition as data frames.

# The next loop will iterate over the programming seeds represented by
# integers which will determine a random yet reproducible set of sample indices
# for each seed environment. The distinct subset of samples will imply
# different feature importance results when the data set size is small, as
# it occurs in this case. An average importance result can be computed in
# order to obtain the most realistic and robust conclusions.

for(i in 1:seeds){
  
  set.seed(i) # A programming seed is established.
  ipartitionIndices <- as.numeric(createDataPartition(y = as.factor(noISclassesVector),
                                                      p = partitionProportion,
                                                      list = FALSE))
  # generates a set of random sample indices determined by the set seed.
  
  icholestasisTrainingData <- t(medianNorm[, ipartitionIndices]) # separates
  # the chosen subset of samples from the full data set and now samples are
  # disposed by rows and columns represent features.
  icholestasisTrainClasses <- noISclassesVector[ipartitionIndices] # subset of
  # class tags that correspond to ordered training samples classes.
  allPartitionedSamples[[paste("seed.", as.character(i), sep = "")]][["training"]] <-
    ipartitionIndices # saves training samples indices from this partition.
  allPartitionedSamples[[paste("seed.", as.character(i), sep = "")]][["test"]] <-
    setdiff(x = 1:length(noISclassesVector), y = ipartitionIndices)
  # saves test samples indices from this partition.
  
  cat(paste("\nCurrent seed: ", as.character(i), ".\n", sep = ""))
  
  set.seed(i) # Again the same programming seed is established just before
  # executing the random forest importance-based feature selection algorithm.
  ifortunaCholestasis <- fortuna(dataset = icholestasisTrainingData,
                                 groups = as.factor(icholestasisTrainClasses),
                                 treeNumber = 1000,
                                 retentionByProp = FALSE,
                                 retentionProportion = 0.8,
                                 retentionLoops = 30,
                                 chosenIterationsNumber = 30,
                                 isBinaryApproach = FALSE,
                                 chosenFeatsClassMethod = "k", # "k", "elbow", "p"
                                 chosenPValue = 0.05,
                                 pAdjMethod = "fdr",
                                 kFeats = nFirstMostImportant,
                                 isElbowPlot = TRUE)
  # the user-defined random forest importance-based feature selection algorithm
  # is executed for each data partition. For each results set we will only
  # consider the most important proteins from which find averaged robust
  # conclusions.
  
  allKBestProts[[paste("seed.", as.character(i), sep = "")]] <-
    ifortunaCholestasis[["selected.features"]] # the most important
  # proteins are recorded and saved for each iteration.
  iorder <- rownames(ifortunaCholestasis[["importance.matrix"]][["impArr"]])
  # All proteins order is retrieved.
  dataFramesList[[paste("seed.", as.character(i), sep = "")]] <-
    cbind.data.frame(proteomeDescription[iorder,],
                     ifortunaCholestasis[["importance.matrix"]][["impArr"]])
  # A data frame considering proteins descriptions and importance results is
  # built.
  
  save(x = allPartitionedSamples, file = "rda/allPartitionedSamples.rda")
  # A record of the partitions sample indices is saved as an RDA file.
  save(x = allKBestProts, file = "rda/allKBestProts.rda")
  # A record of the most important proteins for each partition is saved as an
  # RDA file.
  save(x = dataFramesList, file = "rda/dataFramesList.rda")
  # A record of the feature importance results for each partition is saved as
  # an RDA file.
  
}

load(file = "rda/allPartitionedSamples.rda") # loads entire sample partitioning
# record.
load(file = "rda/allKBestProts.rda") # loads most important features by
# partition record.
load(file = "rda/dataFramesList.rda") # loads entire feature importance results
# tables record.

protRecurrenceArr <- matrix(data = NA, # All data will be empty values at start.
                            nrow = length(allKBestProts), # The number
                            # of rows will equal the number of all partitions.
                            ncol = length(unique(unlist(allKBestProts))),
                            # The number of columns will be equal to the number
                            # of the proteins contained in the union of all the
                            # most important proteins from all data partitions.
                            dimnames = list(names(allKBestProts),
                                            unique(unlist(allKBestProts)))
                            # Rows will be named after a partition index, and
                            # columns will be named using proteins IDs.
                            )
# This variable represents a shell that will be completed based on the
# recurrence of selected features as the most important ones from each
# partition across all partitions. In a matrix representing each partition as a
# row and each column as a feature that at least was selected once among the
# most important ones in any partition, a feature selection recurrence can be
# traced using for instance a binomial test to highlight those proteins that
# seem to be more important regardless specific partitions.

# The number of columns in the previous matrix equals 307 which, after 100
# different partitions from a data set of 40 samples in total, reveals that
# conclusions derived from a subset of samples could induced a skewed view
# of what is the biological events that are most biologically meaningful;
# averaging these results will secure robust deductions.

for(i in 1:nrow(protRecurrenceArr)){ # iterates over all row indices from the
  # previous matrix shell for protein selection recurrence analysis.
  
  protRecurrenceArr[i,] <- as.numeric(colnames(protRecurrenceArr) %in% allKBestProts[[i]])
  # computes a series of binary values that corresponds to the proteins
  # coinciding the matrix column names that have been chosen as most important
  # ones for any specific partition or, by contrasts, were not chosen as so.
  # The more times a protein was considered important, the more statistically
  # significant will be.
  
}

# Last protein selection matrix concept can be depicted as a colored square
# board:

par(mar = c(5.1, 4.1, 4.1, 2.1)) # adjusts R plot frame.

plot(x = protRecurrenceArr, # Matrix object to be plotted.
     col = c("royalblue3", "red"), # Chosen colors to represent matrix values.
     na.col = "black", # In case there are values left, they will be depicted
     # as black squares.
     main = paste("Matrix of selected proteins as hits as part of the ",
                  as.character(nFirstMostImportant),
                  " most important proteins for each partition",
                  sep = ""), # Image title.
     xlab = "Protein", # Horizontal axis name label.
     ylab = "", # Vertical axis name label.
     axis.row = NULL, # refers to rows axis labels size and tilting.
     axis.col = list(cex.axis = 0.55, las = 2),
     # refers to columns axis labels size and tilting.
     key = NULL # refers to the image legend.
     )

# With the next code line, a binomial p-value is applied based on feature
# selection after all data partitions.

ps <- apply(X = protRecurrenceArr, # Calculations will be applied over the
            # elements of a matrix.
            MARGIN = 2, # Those elements correspond to the matrix columns.
            FUN = function(x) binom.test(x = sum(x), # Number of times a
                                         # protein has been selected across all
                                         # partitions.
                                         n = length(x), # Constant total number
                                         # of partitions.
                                         p = sum(protRecurrenceArr) / prod(dim(protRecurrenceArr)),
                                         # represents binomial test hit constant
                                         # probability.
                                         alternative = "greater" # Only
                                         # features that have been selected a
                                         # significantly high number of times
                                         # will appear as significant, and
                                         # those features that have not been
                                         # selected will not, so the probability
                                         # density function will only consider
                                         # significantly high values.
                                         )$p.value)

ps <- ps[order(ps, decreasing = FALSE)] # The binomial p-value vector is
# ordered from most significant features to least.
adj_ps <- p.adjust(p = ps, method = "fdr") # adjusts binomial test protein
# p-values.

write.xlsx(x = cbind.data.frame(proteomeDescription[names(ps),],
                                "binom.pvalue" = ps,
                                "adj.pvalue" = adj_ps),
           # Assembled data frame containing proteins names, description, gene
           # name and binomial test p-value results, with their associated
           # p-values, both unadjusted and adjusted.
           file = "fortuna/recurrencePValues.xlsx", # path of the file to be
           # created.
           col.names = TRUE, # Column names will appear in the file.
           row.names = FALSE, # Row names will not appear in the file.
           sheetName = "Sheet1", # Name of the sheet of the file containing
           # output data.
           append = FALSE # A new file will be generated. A file sharing the
           # same path will be overwritten.
           )

thres <- 0.05 # statistical significance threshold.
top20AdjProts <-
  as.character(na.omit(names(adj_ps[adj_ps <= thres])[1:nFirstMostImportant]))
# retrieves the 20 most important significant proteins from previous analysis.
# This set of proteins can be used to design machine learning classifiers
# that can show a high classification efficacy to group patients from the
# liver diseases that these patients samples pose.

save(x = ps, file = "rda/ps.rda") # saves importance feature selection binomial
# test p-values as an RDA file.
save(x = adj_ps, file = "rda/adj_ps.rda") # saves importance feature selection
# binomial test adjusted p-values as an RDA file.
save(x = top20AdjProts, file = "rda/top20AdjProts.rda") # saves importance
# feature selection binomial test 20 most significant features as an RDA file.

# To test selected significant features classification power, a new data
# partition will be run (it will be reproducible by the seed number 101).

load(file = "rda/medianNorm.rda") # loads median-normalized data.
load(file = "rda/noISclassesVector.rda") # loads classes tags vector.
load(file = "rda/top20AdjProts.rda") # loads 20 most significant proteins after
# importance feature selection binomial test results.

classColName <- "patient" # classes column name for data frame that will
# contain significant features data and samples classes tags.
topFormula <- as.formula(paste(classColName,
                               " ~ ",
                               paste(top20AdjProts, collapse = " + "),
                               sep = "")) # formula object that relates class
# feature to significant features.

partitionProportion <- 0.7 # Proportion of data to be separated as training set.
set.seed(seeds + 1) # establishes a new seed environment.
seed101Indices <- as.numeric(caret::createDataPartition(y = as.factor(noISclassesVector),
                                                        p = partitionProportion,
                                                        list = FALSE))
# A new data set partition is produced to generate new samples indices to build
# a different training data subset.

topTrainingClasses <- as.factor(noISclassesVector)[seed101Indices] # retrieves
# training subset class tags.
topTrainingData <- cbind.data.frame(t(medianNorm[top20AdjProts, seed101Indices]),
                                    "class" = topTrainingClasses)
# builds data frame holding the significant features data and the samples class
# feature tags.
colnames(topTrainingData)[ncol(topTrainingData)] <- classColName # incorporates
# class tags data column name.
topTestClasses <- as.factor(noISclassesVector)[-seed101Indices] # does the same
# class feature tags retrieval but for the test samples.
topTestData <- t(medianNorm[top20AdjProts, -seed101Indices]) # separates
# test subset samples data of the most important features.

# From these data subsets, a collection of machine learning models will be
# built. These models will be different from each other by the classification
# algorithm their based upon. The cross-validation used method will be
# constant, which will be leave-one-out cross-validation.

allClassifyingAlgorithms <- c("rf", "xgbTree", "lda", "knn",
                              "nb", "svmRadial", "glmnet") # caret's terms used
# to refer to the algorithms that want to be executed to produce different
# classification models.
completeAlgorithmNames <- c("Random Forest",
                            "Extreme Gradient Boosting",
                            "Linear Discriminant Analysis",
                            "K-Nearest Neighbors",
                            "Naive Bayes",
                            "Support Vector Machines",
                            "Logistic Regression") # A character vector with
# the complete name of the machine learning algorithms used in this analysis.
# The algorithms order corresponds to the same order the previous vector with
# caret's machine learning algorithm terms has.

topClassifiers <- list() # An empty list that will hold produced classifiers
# results.

# A loop will be run to collect all the results coming from the different
# algorithms. These results will include training samples classification
# performance, test set classification performance and those data needed to
# plot images that represent such classification performance (for instance, 
# ROC curve figure data).

for(i in 1:length(allClassifyingAlgorithms)){ # This loop iterates over all
  # of the caret's machine algorithm character vector elements.
  
  ialgorithm <- allClassifyingAlgorithms[i] # selects one of the algorithms.
  
  if(ialgorithm == "rf"){ # If the machine learning classifier expression
    # matches to one that represents random forest:
    set.seed(1) # establishes reproducible code results. Any time a
    # cross-validation method is applied, results with caret can differ, so
    # a programming seed is set to make sure computational events can be
    # exactly replicated.
    itopClassifier <- train(topFormula, # Formula indicating the class feature
                            # and the forecasting features.
                            data = topTrainingData, # Data frame containing
                            # class feature and significant prediction proteins
                            # data.
                            method = ialgorithm, # Classification algorithm
                            # that is going to be applied.
                            trControl = trainControl(method = "LOOCV",
                                                     classProbs = TRUE),
                            # Constant cross-validation method (leave-one-out
                            # cross-validation, LOOCV) is chosen.
                            ntree = 1000 # argument that is specifically used
                            # when training a random forest with caret. It
                            # represents the number of decision trees a random
                            # forest is made of.
                            )
  }
  if(ialgorithm == "glmnet"){ # If the machine learning classifier expression
    # matches to one that represents logistic regression (elastic net):
    set.seed(1) # A programming seed is set.
    itopClassifier <- train(topFormula, # Formula object so the prediction is
                            # based on the significant features.
                            data = topTrainingData, # Data frame with all
                            # training samples data.
                            method = ialgorithm, # Logistic regression
                            # classification algorithm will be applied.
                            trControl = trainControl(method = "LOOCV",
                                                     classProbs = TRUE),
                            # Constant cross-validation method (leave-one-out
                            # cross-validation, LOOCV) is chosen.
                            family = "multinomial" # This argument refers to
                            # the specific logistic regression mathematical
                            # formula the algorithm will be based upon. It is
                            # the only possible argument value to be used due
                            # to the multiclass classification nature of the
                            # analysis.
                            )
    # A logistic regression model using leave-one-out cross-validation is
    # built.
  }
  if(ialgorithm != "rf" & ialgorithm != "glmnet"){ # If the machine learning
    # classifier expression matches one that represents a different
    # algorithm from random forest or logistic regression:
    set.seed(1) # A programming seed is set for pseudorandom reproducible
    # environment.
    itopClassifier <- train(topFormula,
                            data = topTrainingData,
                            method = ialgorithm,
                            trControl = trainControl(method = "LOOCV",
                                                     classProbs = TRUE))
    # any other classification algorithm is executed to train a model.
  }
  
  itopTestPredClass <- predict(object = itopClassifier, # Classification model.
                               newdata = topTestData, # Test subset data.
                               type = "raw" # Class prediction will be computed.
                               )
  # Using previously trained model, classes for the test data subset are
  # predicted.
  
  itopTestPredProbs <- predict(object = itopClassifier, # Classification model.
                               newdata = topTestData, # Test subset data.
                               type = "prob" # Probabilities for each class
                               # belonging are computed
                               )
  # Using previously trained model, probabilities for each class belonging for
  # the test data subset are predicted.
  
  itopConfRes <- confusionMatrix(data = itopTestPredClass, # Test subset
                                 # samples predicted classes.
                                 topTestClasses # Test subset samples real
                                 # classes. 
                                 )
  # Both test set real classes and test set predicted classes are compared
  # between each other in a confusion matrix. This object also has information
  # about other classification efficacy measurements such as class sensibility,
  # class specificity, balanced accuracy...
  
  # ROC curves are one of the most popular representations of classification
  # efficacy measurement in the machine learning field. "verification" package
  # allows ROC curve plotting based on easy-to-fill arguments to generate a
  # prorated ROC curve. ROC curves have been plotted following the multiclass
  # "one vs. rest" strategy.
  
  topObservedROC <- c() # Empty vector that will hold binary values
  # representing samples class belonging in order.
  topPredROC <- c() # Empty vector that will hold probability values
  # representing samples class belonging in order.
  
  for(j in 1:nlevels(topTrainingClasses)){
    
    jclass <- as.character(levels(topTrainingClasses)[j]) # A class name is
    # selected from all the unique classes in this analysis.
    jobserved <- as.numeric(as.character(topTestClasses) == jclass)
    # Binary vector generated from a boolean vector that represents if the
    # currently analysed class corresponds to the class tags vector elements.
    topObservedROC <- append(x = topObservedROC, values = jobserved)
    # appends binary vector to elongate classes cumulative vector.
    jpredicted <- itopTestPredProbs[,jclass] # retrieves test set samples
    # probabilities of having been identified as the current analysed class.
    topPredROC <- append(x = topPredROC, values = jpredicted) # appends
    # probability values to cumulative prediction probability vector.
    
  }
  
  itopROC <- verification::roc.plot(x = topObservedROC, # Elongated binary
                                    # vector referring to which class tags
                                    # coincide with vector class tags by
                                    # indicating hits with 1 and non-hits
                                    # with 0.
                                    pred = topPredROC, # A prediction
                                    # probability assign to each sample for
                                    # each class.
                                    xlab =  "1 - specificity", # Horizontal
                                    # axis label.
                                    ylab = "Sensitivity", # Vertical axis
                                    # label.
                                    show.thres = FALSE # shows critical
                                    # discerning ROC curve thresholds if it
                                    # equals "TRUE".
                                    )
  # generates ROC curve plot and its related information. It even contains
  # all of its graphical parameters to reproduce it on a plot and they allow
  # you to change its appearance.
  
  topClassifiers[[ialgorithm]] <- list("classifier" = itopClassifier,
                                       "testPred" = itopTestPredClass,
                                       "testProb" = itopTestPredProbs,
                                       "confArr" = itopConfRes,
                                       "ROC" = itopROC) # collects all
  # classification performance parameters into a single list that will be
  # associated to the used classification algorithm.
  
  save(x = topClassifiers, file = "rda/topClassifiers.rda") # saves and updates
  # classifiers results.
  
}

load(file = "rda/topClassifiers.rda") # loads results from all classifiers
# based on most significant proteins.

topAUCs <- c() # Empty vector to hold the set of area under the curve values
# from all ROC curves to eventually plot them as part of a legend.
ROClegend <- c() # Variable which will hold a character vector to provide
# a legend for a multiple ROC curve image.
paletteForAlgorithms <-
  RColorBrewer::brewer.pal(n = 8, name = "Set2")[1:length(allClassifyingAlgorithms)]
# selects a combination of colors to be used in the ROC curves plot.
names(paletteForAlgorithms) <- allClassifyingAlgorithms
# names previous color vector to associate each color to a classification
# algorithm.

bestAlgConfArr <- names(which.max(unlist(lapply(X = topClassifiers,
                                                FUN = function(x)
                                                  x$ROC$roc.vol$Area))))
# retrieving best classifier name based on the area under the ROC curve
# parameter according to test data set classification efficacy.

for(i in 1:length(topClassifiers)){ # This loop will iterate over all ROC curves
  # associated information considering each model results to plot their image.
  
  ialgorithm <- names(topClassifiers)[i] # selects the algorithm that was used
  # to generate some specific ROC curve results.
  iroc <- topClassifiers[[ialgorithm]][["ROC"]] # accesses to its attached ROC
  # curve information.
  
  xROC <- rev(matrix(iroc[["plot.data"]], ncol = ncol(iroc[["plot.data"]]))[, 3])
  # retrieves horizontal coordinate information to plot the ROC curve and
  # reverses it.
  yROC <- rev(matrix(iroc[["plot.data"]], ncol = ncol(iroc[["plot.data"]]))[, 2])
  # retrieves vertical coordinate information to plot the ROC curve and
  # reverses it.
  newXROC <- c() # Empty vector that will hold transformed horizontal axis plot
  # values. The reason why this part was implemented is that some algorithms,
  # due to the class prediction probability they generate (e. g., KNN),
  # sometimes they can produce a ROC curve without a staircase appearance. This
  # code chunk just makes sure ROC curves have that form without altering
  # their area under the curve value.
  newYROC <- c() # Empty vector that will hold transformed vertical axis plot
  # values. This adaptation allows a ROC curve staircase representation,
  # considering that some classifiers produce a set of class belonging
  # probabilities that alters ROC curve appearance.
  
  for(j in 1:length(xROC)){ # Loop that iterates over all ROC horizontal
    # coordinates values.
    
    jx1 <- xROC[j] # orderly isolates one horizontal axis coordinates value.
    jy1 <- yROC[j] # does the same thing as previous code line, but over the
    # vertical axis coordinates. The length of both coordinates vectors is the
    # same which is mandatory to eventually generate a plot.
    
    if(j != 1){ # If that loop has got over its first element, this part of the
      # code will be executed.
      
      jx0 <- xROC[j - 1] # Previous horizontal axis coordinates vector element
      # is saved.
      jy0 <- yROC[j - 1] # Previous vertical axis coordinates vector element is
      # saved.
      jxBool <- jx0 == jx1 # Current vector element and previous vector element
      # both coming from the horizontal axis vector are compared to determine
      # if they are equal or not.
      jyBool <- jy0 == jy1 # Current vector element and previous vector element
      # both coming from the vertical axis vector are compared to determine if
      # they are equal or not.
      
      if(as.numeric(jxBool) + as.numeric(jyBool) != 1){
        # As long as both are "FALSE" or both are "TRUE", this will be executed.
        
        newXROC <- append(x = newXROC, values = mean(c(xROC[j - 1], xROC[j])))
        # The mean value of previous and current horizontal axis coordinates
        # elements is calculated and appended to transformed horizontal axis
        # values vector.
        newYROC <- append(x = newYROC, values = yROC[j - 1])
        # Previous vertical axis coordinates vector element is incorporated to
        # transformed vertical axis values vector.
        
        newXROC <- append(x = newXROC, values = mean(c(xROC[j - 1], xROC[j])))
        # The mean value of previous and current horizontal axis coordinates
        # elements is calculated and appended to transformed horizontal values
        # vector.
        newYROC <- append(x = newYROC, values = yROC[j])
        # Current vertical axis coordinates vector element is incorporated to
        # transformed vertical axis values vector.
        
      }
      
    }
    
    newXROC <- append(x = newXROC, values = xROC[j]) # Next, current
    # horizontal axis coordinates vector element is added to transformed
    # horizontal axis values vector.
    newYROC <- append(x = newYROC, values = yROC[j]) # Current vertical axis
    # coordinates vector elements are added to transformed vertical axis values
    # vector.
    
  }
  
  if(ialgorithm == bestAlgConfArr){ # If we are handling the set of values
    # coming from the best working classifier during a loop iteration:
    newXROC_best <- newXROC # ROC curve horizontal coordinates values are saved
    # to represent the ROC curve later.
    newYROC_best <- newYROC # ROC curve vertical coordinates values are saved
    # to represent the ROC curve later.
  }
  
  if(i == 1){ # If the loop is executing its first iteration, "plot" function
    # needs to be run before "points" function, used to create images over
    # an already plotted graph.
    
    plot(x = newXROC, # Transformed ROC curve horizontal axis coordinates
         # values.
         y = newYROC, # Transformed ROC curve vertical axis coordinates
         # values.
         xlim = range(pretty(newXROC)), # ranges horizontal values to try
         # to represent aesthetic horizontal axis limits.
         ylim = range(pretty(newYROC)), # ranges vertical values to try
         # to represent aesthetic vertical axis limits.
         yaxt = "n", # Vertical axis numerical labels will not be plotted.
         type = "l", # A linear type of plot will be represented.
         col = paletteForAlgorithms[i], # will color the plotted ROC curve.
         lwd = 3, # Plot line width.
         main = "Algorithms performance ROC curves", # Plot title.
         xlab = "1 - specificity", # Horizontal axis title.
         ylab = "sensitivity", # Vertical axis title.
         cex.axis = 1.5, # Size of numerical labels from horizontal and
         # vertical axes.
         cex.lab = 1.5, # Size of axes titles labels.
         cex.main = 2 # Size of plot main title.
         ) # plots ROC curve image.
    axis(side = 2, # refers to conventional vertical left axis.
         las = 2, # controls axis numerical tags orientation. "2" corresponds
         # to readable text orientation.
         at = pretty(newYROC), # generates a set of linearly spaced values
         # to make them those representing the vertical axis marks.
         labels = pretty(newYROC), # generates a set of linearly spaced values
         # to make them those representing the vertical axis labels.
         cex.axis = 1.5 # controls axis numbers size.
         ) # plots vertical axis number tags.
    
  }else{ # If there is already a plot offering a frame to represent the rest of
    # ROC curves, this part of the code will be executed.
    
    points(x = newXROC, # Transformed ROC curve horizontal axis coordinates
           # values.
           y = newYROC, # Transformed ROC curve vertical axis coordinates
           # values.
           type = "l", # Linear type of plot.
           col = paletteForAlgorithms[i], # Chosen color for ROC curve.
           lwd = 3 # Line width.
           )
    
  }
  
  topAUCs <- append(x = topAUCs, values = iroc$roc.vol$Area) # compiles the
  # area under the ROC curve values from each classification model.
  ROClegend <- append(x = ROClegend,
                      values = paste("AUC (",
                                     completeAlgorithmNames[i],
                                     "): ",
                                     as.character(round(topAUCs[i], 2)),
                                     sep = "") # adds a legend line showing a
                      # classification algorithm and its related area under
                      # the ROC curve result.
                      )
  
}

names(topAUCs) <- allClassifyingAlgorithms # names each area under the ROC
# curve vector value after the classification algorithm that produced that
# value.
names(ROClegend) <- allClassifyingAlgorithms # names each character vector
# element relating a classification algorithm to an area under the curve result
# using the name of the caret term utilized to execute each classification
# algorithm.

legend(x = "bottomright", # Legend position in the plot.
       inset = c(0.01, 0.01), # Legend position placement correction.
       legend = ROClegend, # Legend text.
       pch = 22, # Symbols that are wanted to be used in the legend.
       col = "black", # Color of symbols boxes.
       pt.bg = paletteForAlgorithms, # Color of symbols insides.
       cex = 1.5, # controls legend size.
       pt.cex = 2 # controls legend symbols size.
       ) # plots legend for the ROC curves image.

plot(x = newXROC_best, # Horizontal axis ROC curve transformed values
     # representing the best classifier based on its area under the ROC curve
     # result.
     y = newYROC_best, # Vertical axis ROC curve transformed values
     # representing the best classifier based on its area under the ROC curve
     # result.
     xlim = range(pretty(newXROC_best)), # establishes horizontal axis upper
     # and lower limit.
     ylim = range(pretty(newYROC_best)), # establishes vertical axis upper
     # and lower limit.
     yaxt = "n", # Vertical axis values will not be directly plotted through
     # this function.
     type = "l", # makes the plot one which coordinates are chained as a line.
     col = paletteForAlgorithms[bestAlgConfArr], # chooses a color based on
     # the chosen algorithm model which performed best to plot its related
     # ROC curve.
     lwd = 3, # determines line width.
     main = "Performance ROC curve", # includes a plot title.
     xlab = "1 - specificity", # Horizontal axis name label.
     ylab = "sensitivity", # Vertical axis name label.
     cex.axis = 1.5, # controls numerical axis labels size.
     cex.lab = 1.5, # controls axes names labels size.
     cex.main = 2 # determines plot title size.
     )
# plots ROC curve with the greatest area under the curve associated value.
axis(side = 2, # refers to conventional vertical left axis.
     las = 2, # controls axis numerical tags orientation. "2" corresponds
     # to readable text orientation.
     at = pretty(newYROC_best), # generates a set of linearly spaced values
     # to make them those representing the vertical axis marks under which its
     # numerical tags will be placed.
     labels = pretty(newYROC_best), # generates a set of linearly spaced values
     # to make them those representing the vertical axis labels.
     cex.axis = 1.5 # controls axis numbers size.
     ) # plots vertical axis numbers labels.

legend(x = "bottomright", # Legend position.
       inset = c(0.01, 0.01), # Legend position correction (horizontally and
       # vertically).
       legend = ROClegend[bestAlgConfArr], # retrieves name of the
       # classification algorithm and area under the ROC curve results from
       # character vector compiling all classification algorithms results.
       pch = 22, # Used legend symbol type.
       col = "black", # Legend symbols frame line color.
       pt.bg = paletteForAlgorithms[bestAlgConfArr], # Legend symbols interior
       # color, which will match the ROC curve line one.
       cex = 2, # controls legend size.
       pt.cex = 2.5 # controls legend symbols size.
       ) # incorporates legend for best performing classifier ROC curve plot.

testTableResults <- matrix(data = unlist(lapply(X = topClassifiers,
                                                FUN = function(x) x[["confArr"]][["overall"]][c("Accuracy", "Kappa", "AccuracyPValue")])),
                           ncol = length(allClassifyingAlgorithms))
# Using a function from the so-called "apply" family, considering the list
# structure of the object compiling all classifier results over the test subset
# samples, it retrieves the accuracy, kappa coefficient and accuracy p-value
# associated to each algorithm and disposes them as a table.
testTableResults <- rbind(testTableResults, topAUCs) # orderly incorporates
# the area under the curve values to the test set results table.
testTableResults <- round(testTableResults, 4) # rounds all numbers to make
# them show only 4 decimal places.
rownames(testTableResults) <- c("accuracy", "kappa", "acc.pvalue", "AUC")
# names test set results table rows after the performance parameters they
# represent.
colnames(testTableResults) <- c("RF", "XGB", "LDA", "KNN", "NB", "SVM", "LR")
# names test set results table rows after the classification algorithms that
# were used.

write.xlsx(x = testTableResults, # Data table to be converted into Excel file.
           file = "classifierResults.xlsx", # Chosen file name.
           row.names = TRUE, # Table row names will be added to file.
           col.names = TRUE, # Table column names will be added to file.
           append = FALSE, # Information will overwrite files with the same
           # name.
           sheetName = "Sheet1" # Excel file sheet to be written.
           ) # writes test set classification performance parameters in the
# form of Excel file.

load(file = "rda/proteomeDescription.rda") # loads RDA file containing protein
# descriptions table.
load(file = "rda/medianNorm.rda") # loads RDA file containing normalized
# protein data abundance table.

medianNorm_gene <- t(medianNorm) # transposes normalized protein abundance data
# matrix. A transposed copy is being created to overwrite some column names
# and row names of the matrix to facilitate some preferred protein nomenclature
# and sample representation from this array.
colnames(medianNorm_gene) <- proteomeDescription[colnames(medianNorm_gene), "gene.ID"]
# names proteins using their corresponding gene name.

PFICtagsFile <- "originalFiles/Tags_PFIC.xlsx" # path and name of an additional
# file containing the information of the subtypes of PFIC the patients from
# which their samples come.

# It will be interesting to plot a heat map considering those significant
# proteins in the classification task and see how different PFIC subtypes
# are grouped in the image.

PFICtags <- read.xlsx2(file = PFICtagsFile, sheetIndex = 1) # loads Excel file
# containing patients' PFIC subtypes table.

# Using this information, the normalized protein data matrix will be updated to
# consider also PFIC subtypes by overwriting samples names. A number
# showing the PFIC subtype will be included at the end of the samples names.

for(i in 1:nrow(PFICtags)){ # This loop iterates over the loaded table
  # having information about the PFIC subtypes of patients' samples.
  
  oldPFICsampleName <- paste(PFICtags[i, "TMT"],
                             PFICtags[i, "Tag"],
                             "PFIC",
                             sep = ".") # generates one of the PFIC samples
  # names starting from the PFIC subtypes table information. 
  
  boolIndex <- rownames(medianNorm_gene) == oldPFICsampleName
  # locates previous sample name in normalized protein abundance data matrix
  # row names. It has to be noted that this matrix we are working with right
  # now corresponds to that that was produced as a copy.
  
  newPFICsampleName <- paste(PFICtags[i, "TMT"],
                             PFICtags[i, "Tag"],
                             PFICtags[i, "PFIC.subtype"],
                             sep = ".") # assembles new sample name having
  # the PFIC subtype number.
  
  rownames(medianNorm_gene)[boolIndex] <- newPFICsampleName
  # substitutes old PFIC sample name by its new name with the PFIC subtype
  # indicated.
  
}

medianNorm_gene <- medianNorm_gene[c(1:28,
                                     37,38, # Samples subtype: PFIC1.
                                     34:36, # Samples subtype: PFIC2.
                                     29:33, # Samples subtype: PFIC3.
                                     39:40), # Samples subtype: PFIC4.
                                   ] # manually disposes PFIC samples by their
# subtype in a wanted order.
dataForPlots <- medianNorm_gene[, proteomeDescription[top20AdjProts,"gene.ID"]]
# indexes matrix to just select those proteins involved in correct machine
# learning classification of samples that were detected as significant.

namedColorVector <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(3,1,2,5)],
                      "gold") # A color vector containing the same colors used
# to depict the classes before is generated, and in this case it will also
# be named, so each color will have the name of its class attached to it.
names(namedColorVector) <- c("Control", "AATD", "Alagille", "Atresia", "PFIC")
# adds classes as names to classes colors.

set.seed(1) # sets a programming seed. Since a calculation depending on
# randomness (which is random forest) is about to be executed, the set seed
# will help to reproduce the exact results anytime.
rfForPlot <- randomForest(x = dataForPlots, # Subset of significant proteins.
                          # All samples are contained in this table.
                          y = as.factor(noISclassesVector), # Classes vector
                          # according to samples ordering in the significant
                          # proteins data table.
                          proximity = TRUE, # Proximity measurement will be
                          # emitted as part of the result of the classifier.
                          # This allows multidimensional scaling (MDS) plot
                          # representation.
                          importance = TRUE,# Importance measurement for
                          # feature classification power study is recorded.
                          ntree = 1000 # Number of trees the random forest
                          # will have.
                          ) # A classifier using the specialized random forest
# R library is generated. This classifier is useful to make some
# representations that allow understanding the characteristics of the classes
# that are being classified.
RFClassEx(rfObject = rfForPlot, # Random forest object.
          trueClassesOrder = c("Control", "AATD", "Alagille", "Atresia", "PFIC"),
          # Vector that establishes an order in bar parts representation.
          classesColors = namedColorVector[c(c("Control", "AATD", "Alagille", "Atresia", "PFIC"))],
          # Colors that will be matched by their order to each one of the case
          # classes.
          cexMain = 2, # Title size factor.
          cexAxisNumbers = 1.3, # Axes numerical labels size factor.
          cexPlotAxis = 1.2, # Axes name labels size factor.
          tiltAngle = 40, # Axes labels marks orientation factor (it determines
          # the label inclination degrees).
          cexLabels = 1.2, # Features names labels size factor.
          cexLeg = 2, # Figure legend size factor.
          legInset = c(0.02, 0.02), # Figure legend position correction factor.
          adjText = c(1, 1) # Features names labels position correction factor.
          ) # User-defined function that represents feature importance as a bar
# chart to see how valuable each feature is for random forest classification
# efficacy and also how each feature is involved to classify samples from
# specific classes particularly.

reducedRowNames <- unlist(lapply(X = strsplit(x = rownames(dataForPlots),
                                              split = ".",
                                              fixed = TRUE),
                                 FUN = function(x) x[length(x)]))
# splits sample names into multiple chunks that separately indicate their
# preparation batch, Tandem Mass Tag labels and disease class. Subsequently,
# only the class information is extracted. This will be used to name array
# columns to plot the array as a heat map.

dataForPlots_copy <- t(dataForPlots) # transposes significant proteins matrix
# and saves it as a variable.
colnames(dataForPlots_copy) <- reducedRowNames # overwrites column names
# which were the samples names with only the name of the class they can be
# grouped into.

heat(dataF = dataForPlots_copy, # Data matrix containing significant features
     # data and renamed samples using only their class belonging.
     distance = "correlation", # Distance measurement using in heat map
     # dendrograms.
     xcex = 1.3, # Horizontal axis figure labels size factor.
     ycex = 1.2, # Vertical axis figure labels size factor.
     ylab = TRUE, # Boolean factor that determines if vertical axis labels are
     # wanted in the figure.
     rowD = TRUE, # Boolean factor that indicates if distance dendrogram
     # hierarchical clustering among data array rows needs to be computed.
     colD = FALSE, # Boolean factor that indicates if distance dendrogram
     # hierarchical clustering among data array columns needs to be computed.
     colSet = "greenred", # Color palette that has been chosen for the 
     # heat map.
     rotateXLabs = 90, # Degrees the horizontal axis labels will be
     # rotated.
     rotateYLabs = NULL, # Degrees the vertical axis labels will be
     # rotated.
     ADJCOL = c(1, 0.4), # Factors that correct horizontal axis labels position.
     ADJROW = c(0.2, 0) # Factors that correct vertical axis labels position.
     ) # user-defined function that plots a heat map using a function from
# gplots library to make heat maps.

# The next code chunk is the same as the previous code lines, which were used
# to plot a heat map of the significant proteins in the analysis, but in this
# case, a single argument has changed: The next code lines will produce a heat
# map and now a hierarchical clustering dendrogram will also be produced by
# grouping the data considering the columns of the data matrix.

heat(dataF = dataForPlots_copy,
     distance = "correlation",
     xcex = 1.3,
     ycex = 1.2,
     ylab = TRUE,
     rowD = TRUE,
     colD = TRUE,
     colSet = "greenred",
     rotateXLabs = 90,
     rotateYLabs = NULL,
     ADJCOL = c(1, 0.4),
     ADJROW = c(0.2, 0))

# One of the best classifiers by their test subset classification performance
# was random forest. Random forest allows sample similarity calculation by
# comparing the classification proximity each sample has with respect to the
# rest by knowing the leaf node each sample fell into across decision trees.
# Recurring leaf nodes for multiple samples would imply higher similarity
# between them.

distance <- cmdscale(1 - rfForPlot$proximity, eig = TRUE, k = 3) # It executes
# a calculation known as multidimensional scaling (MDS). By subtracting random
# forest sample proximity results to one, since all proximity results are
# values between 0 and 1, a distance measure can be obtained. After generating
# a distance data matrix, the MDS analysis can be executed. It will generate as
# many dimensions as indicated by the "k" argument.

distanceData <- distance$points # retrieves the data matrix from the set of
# produced results.
colnames(distanceData) <- c("MDS1", "MDS2", "MDS3") # renames data matrix
# columns.

plot3D(data3D = distanceData, # MDS distance after random forest sample
       # similarity results.
       classVector = noISclassesVector, # All sample classes vector.
       classColors = namedColorVector, # Named character vector of colors after
       # analysis classes.
       main3D = "", # Plot title.
       cex3D = 2, # Point size.
       cexLegend = 4, # Legend size.
       xAxisLimits = c(-0.6, 0.6), # Width axis limits.
       yAxisLimits = c(-0.6, 0.6), # Height axis limits.
       zAxisLimits = c(-0.6, 0.6), # Depth axis limits.
       par3D = 2 # Axis numerical labels and axes name labels size.
       ) # A user-defined function used to plot 3D images. This will represent
# MDS results after random forest classification.

PCAGif() # User defined function that creates a GIF movie by rotating the 3D
# plot around its vertical axis. It saves the GIF file in current working
# directory. This function allows rotation sense to be modified.

dev.new() # creates new plot frame.

par(mfrow = c(4, 5)) # divides plot frame into boxes to plot more images per
# plot frame. The first number indicates the number of subdivision rows and the
# second number refers to the number of subdivision columns.

yLabelsFactor <- 1/13 # A factor that contributes in a calculation to determine
# boxplots class tags vertical position. The greater the value is, the lower
# the tags height will be.
yTitleFactor <- 10/45 # A factor that contributes in a calculation to determine
# title vertical position. The greater the value is, the higher the title
# position will be.

for(i in 1:ncol(dataForPlots)){ # Loop that iterates over each column which
  # represents a protein.
  
  idat <- cbind.data.frame(dataForPlots[,i], as.factor(noISclassesVector))
  # creates a data frame from one protein data and the classes vector to
  # link each sample from a class to its protein data value.
  colnames(idat) <- c(colnames(dataForPlots)[i], "class")
  # renames columns including the protein gene name.
  idat$class <- factor(idat$class, levels = c("Control", "AATD", "Alagille", "Atresia", "PFIC"))
  # Boxplots classes wanted order.
  
  boxplot(idat[,1] ~ idat[,2], # indicates a formula based on the columns
          # of the data frame from which to plot a boxplot. The first column
          # represents the independent values, and the second column represents
          # the groups or classes the values are related to.
          data = idat, # Data frame containing plot information.
          xaxt = "n", # leaves horizontal axis labels blank.
          yaxt = "n", # leaves vertical axis labels blank.
          ylim = range(pretty(idat[,1])), # ranges vertical axis considering
          # its limit values.
          xlab = "", # Horizontal axis name label.
          ylab = "", # Vertical axis name label.
          main = "", # Plot title.
          pch = 21, # Boxplot points type.
          bg = "grey", # Boxplots color.
          outcex = 1.6, # Boxplot points size.
          cex.lab = 1.8 # Boxplot axes names label size.
          )
  axis(side = 1, # plots labels on horizontal axis.
       at = 1:nlevels(idat[,2]), # Axis marks will be plotted on positions from
       # one to five on the horizontal axis.
       labels = rep("", 5) # No labels will be added to the axis.
       ) # It represents labels on a wanted axis and on determined positions.
  # Thanks to this function, first only axis marks will be added. After that,
  # its text will be included to reference the analysis classes.
  axis(side = 2, # plots labels on vertical axis.
       las = 1, # controls labels orientation; labels will be oriented to be
       # read.
       at = pretty(idat[,1]), # Vertical axis numerical labels will be disposed
       # on these positions.
       labels = pretty(idat[,1]), # Vertical axis numerical values.
       cex.axis = 2.2 # Axis labels size.
       ) # plots axis labels.
  
  highest <- range(pretty(idat[,1]))[2] # Vertical axis upper limit.
  lowest <- range(pretty(idat[,1]))[1] # Vertical axis lower limit.
  difference <- highest - lowest # computes difference between upper and lower
  # limit. This helps setting a proportioned position to the incorporated 
  # labels into the boxplot.
  labelsYPosition <- lowest - (yLabelsFactor * difference) # calculates a
  # position for boxplot classes tags.
  titleYPosition <- highest + (yTitleFactor * difference) # calculates a
  # position for boxplot title.
  
  text(x = 1:nlevels(idat[,2]), # Horizontal axis text positions.
       y = labelsYPosition, # vertical axis text positions.
       labels = c("Control", "AATD", "Alagille", "Atresia", "PFIC"),
       # Text labels.
       adj = c(1, 1), # Labels position readjustment coordinates (the first one
       # refers to the horizontal positioning, and the second one to the
       # vertical placement). Values interval is from 0 to 1.
       xpd = NA, # This argument allows plotting text out of the axes-delimited
       # plot space.
       srt = 35, # Number of degrees text will be tilted with respect to the
       # horizontal axis.
       cex = 2.2 # Size of plotted text.
       ) # incorporates class tags.
  
  text(x = mean(1:nlevels(idat[,2])), # finds and sets horizontal axis middle
       # point to plot boxplot title.
       y = titleYPosition, # Vertical axis wanted position for boxplot title.
       labels = colnames(idat)[1], # sets boxplot title as protein gene name.
       adj = c(0.5, 0.5), # Text adjustment factors.
       xpd = NA, # allows plotting text out of the axes-delimited plot space.
       srt = 0, # Number of degrees text will be tilted with respect to the
       # horizontal axis.
       cex = 3.5 # Size of plotted text.
       ) # incorporates vertical axis name label.
  
}

# 6) Extra images output

# One of the best ways to understand what data have to say consists of just
# representing them in a plain but descriptive way. Boxplots usually offer a
# clear image that shows data trends and corroborate intuitions related to data
# patterns across samples and classes, even though they do not stack up
# complete evidence by themselves. 

# A user-defined function to generate a set of boxplots for each feature was
# written to produce thousands of boxplots from entire proteomics data sets.
# This is useful as a first look inside the data or as a partial assessment
# of data suspected tendencies.

# This function, by its particular way of working, requires a matrix of
# statistical significance p-values which contrasts the data from one protein
# between two compared classes. P-values have to be ordered by columns, one
# contrast by column, and contrasts of compared classes follow a specific order
# considering compared classes names. For this analysis, no statistical
# significance assessment will be considered for boxplots images, so a dummy
# matrix full of ones will be used.

# Anyway, statistical significance column names supposed order will be
# preserved. The next code chunk reveals the contrasts column names order
# considering the name of the compared classes expected by the user-defined
# function:

dummyPValues_ncol <- unique(noISclassesVector)[order(unique(noISclassesVector))]
# alphabetically orders analysis classes names. Even though the variable name
# does not match its utility, that is because later this variable content
# will be replaced with other data. For now, it will hold alphabetically
# ordered unique classes.
dummyPValues_colnames <- c() # Empty vector that will hold p-values matrix
# column names.

for(i in 1:length(dummyPValues_ncol)){ # This loops iterates over all classes
  # names.
  
  for(j in (i+1):length(dummyPValues_ncol)){ # This loops iterates over all the
    # classes names except for the one that is being contemplated by the outer
    # loop. This allows creating all possible paired combinations without
    # repetition from the character vector.
    
    if(i != length(dummyPValues_ncol)){ # Only when the outer loop has not
      # reached its final element, these inner loop code lines will keep being
      # executed.
      
      jstr <- paste("adj.pval.",
                    dummyPValues_ncol[j],
                    "-",
                    dummyPValues_ncol[i],
                    sep = "") # ensembles a string from combining two classes
      # names, with the first one being always that which is alphabetically
      # posterior to the one that is at the end of this string.
      dummyPValues_colnames <- append(x = dummyPValues_colnames,
                                      values = jstr) # appends recently
      # constructed string to a vector that will hold all combinations of
      # classes..
      
      # This string is generated to name those columns of a matrix that will
      # hold statistical significance p-values. For other analyses, a function
      # used to compute the results of an ANOVA was defined. This function
      # produces p-value results in a particular way, mainly because of the
      # order R's ANOVA built-in function results are generated. After building
      # this code module, another function that would produce boxplots showing
      # statistical significance between represented groups was written, so
      # the results from one were the input for the other.
      
      # Now what it is being tried to do is to emulate those results from the
      # ANOVA funciton to not produce any statistically significant result to
      # only generate boxplot images, to reuse the boxplots user-defined
      # function. Some possible set of result is being mimicked, just to
      # generate all the boxplots of this data set.
      
    }
    
  }
  
}

dummyPValues_ncol <- length(unique(noISclassesVector)) # Now this variable will
# hold the number of possible paired combinations without repetition of the
# classes from the analysis. Its value will be replaced by another, and for now
# it equals the number of classes that there are in total.
dummyPValues_ncol <- factorial(dummyPValues_ncol) /
  (factorial(dummyPValues_ncol - 2) * 2)
# The number of possible paired combinations without repetition is calculated
# as a result of the formula of the number of "k-element" combinations of "n"
# objects without repetition, where "k" equals 2 and "n" equals 5 (which is the
# total number of classes).

dir.create(path = "allBoxplots") # creates a directory to hold all boxplots
# images to show the protein levels between different groups. This can be
# clearly interpreted and gives the user a tangible sense of the features
# behaviour for the different conditions.

photoBox2(is.already = TRUE, # indicates if it is necessary to create a
          # directory for all the images or if it has been already created. 
          nameDir = NULL, # indicates the path where the directory that will
          # hold the images is located.
          where = "allBoxplots/", # corresponds to the name of the directory
          # that will hold all boxplots images.
          dat = medianNorm, # Data matrix from which the boxplots will be
          # generated. Features are disposed by columns, and samples are
          # disposed across rows.
          chosenQuality = 100, # A value which affects resolution. The
          # "quality" of the JPEG image, as a percentage. Smaller values will
          # give more compression but also more degradation of the image.
          chosenWidth = 1000, # Image width.
          chosenHeight = 800, # Image height. Together with the previous
          # argument, they affect image resolution.
          pvalues = matrix(data = 1,
                           nrow = nrow(medianNorm),
                           ncol = dummyPValues_ncol,
                           dimnames = list(rownames(medianNorm),
                                           dummyPValues_colnames)), # creates
          # a matrix full of ones, representing statistical significance
          # results. Using the argument "dimnames", its columns have been named
          # following the necessary criteria for the function to work.
          groups2 = noISclassesVector, # Complete classes tags vector (as a
          # character vector object).
          is.ANOVA = TRUE, # It indicates if the significance results come
          # from an ANOVA or a t-test. It has become an obsolete argument.
          main3 = reducedDescription(proteomeDescription)[rownames(medianNorm)],
          # generates all boxplot images titles in order.
          isNorm2 = TRUE # Have values been normalized or not?
          ) # creates a set of boxplots images for different groups and one
# image per feature. All images are generated in the specified path. For this
# function to work, it requires a data matrix with its features arranged by
# columns and a p-value array with its columns named in a specific way. Those
# names have to contemplate pairwise class combinations and they have to be
# ordered by their name in a specific fashion.
