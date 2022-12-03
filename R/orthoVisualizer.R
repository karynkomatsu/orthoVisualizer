
# Exported functions================

#' Annotate Ortholog Subsequence (gene/motif) Occurrence
#'
#' A function that prints each sequence in FASTA file and annotates every
#'    occurrence of user-inputted ortholog pattern subsequence.
#'
#' @param fastaPath is a character string indicating the file path for input
#' FASTA file containing DNA sequences.
#' @param pattern is a character string indicating the ortholog gene/motif
#' subsequence to search for .
#'
#' @return Returns invisible NULL. Prints each sequence in the FASTA file with
#'    annotation.
#'
#' @examples
#'
#' #Example 1:
#' #Print every sequence in FASTA file containing 1 DNA sequence and annotate
#' #every occurrences of target orthologous binding motif pattern "CGCG".
#' # Get raw data
#' placenta_path <-  system.file("extdata", "P_Placenta_seq.fasta", package = "orthoVisualizer")
#' annotateSeq(fastaPath = placenta_path, pattern = "CGCG")
#'
#' #Example 2:
#' #Print every sequence in FASTA file containing 6 DNA sequences and annotate
#' #every occurrences of target orthologous binding motif pattern "CGCG".
#' # Get raw data
#' DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
#' annotateSeq(fastaPath = DNA_seq_path, pattern = "CGCG")
#'
#' @export
annotateSeq <- function(fastaPath, pattern) {

  # Read, validate, and load in the FASTA data from path
  fastaDF <- loadFASTA(fastaPath)  # Dataframe for FASTA sequences

  start <- stop <- NULL

  for (i in 1:length(fastaDF$seq)){
    # Get which of  split the upstream into lines of 50 nucleotides
    start <- seq(1, nchar(fastaDF$seq[i]), 50)             # start indices
    stop <- seq(50, nchar(fastaDF$seq[i]), 50)             # stop indices

    # Add the final stop index AKA the final incomplete line
    if (!(nchar(fastaDF$seq[i]) %in% stop)){
      stop <- append(stop, nchar(fastaDF$seq[i]))
    }

    # Add all the lines of sequences
    splitSeq <- base::substring(fastaDF$seq[i], start, stop)

    # Iterate seq and find any motifs that could be the inputted binding motifs
    # (ie. find all matches of the inputted subsequence)
    # Vectors of motifs & starting indices
    motifIndex <- findMotifs(fastaDF$seq[i], pattern)

    # Print formatted genome sequence with corresponding annotations
    printAnnotate(fastaDF$head[i], splitSeq, pattern, motifIndex)
  }
  return(invisible(NULL))
}

#' Count Ortholog (gene/motif) Occurrence
#'
#' Plot Ortholog Subsequence (gene/motif) Occurrence
#'
#' A function that plots how many times user-inputted ortholog pattern appears
#'    in each sequence in FASTA file. Each bar represents one sequence.
#'
#' @param fastaPath is a character string indicating the file path for input
#' FASTA file containing DNA sequences.
#' @param pattern is a character string indicating the ortholog gene/motif
#' subsequence to search for
#'
#' @return Returns a tibble of "assigned id", "header", "number of occurrence",
#'    "length of sequence", "ratio of occurrence (frequency ratio)" in order,
#'    for every sequence in FASTA file
#'
#' @examples
#'
#' #Example 1:
#' #Give the quantity tibble for FASTA file containing 1 DNA sequence (1 bar)
#' #and target orthologous binding motif pattern "CGCG".
#' # Get raw data
#' placenta_path <-  system.file("extdata", "P_Placenta_seq.fasta", package = "orthoVisualizer")
#' quantSeq(fastaPath = placenta_path, pattern = "CGCG")
#'
#' #Example 2:
#' #Give the quantity tibble for FASTA file containing 6 DNA sequences (6 bars)
#' #and target orthologous binding motif pattern "CGCG".
#' # Get raw data
#' DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
#' quantSeq(fastaPath = DNA_seq_path, pattern = "CGCG")
#'
#' @references
#' Muller, K., et al. In-Line Documentation for R
#' R Package Roxygen2 Version 7.2.2. Comprehensive R Archive Network (CRAN),
#' 22 July. 2022, \href{https://cran.r-project.org/web/packages/roxygen2/index.html}{Link}.
#'
#' @export
#' @importFrom tibble tibble
quantSeq <- function(fastaPath, pattern) {

  # Input: Path to FASTA file,
  #        Pattern of the target orthologous motif / gene subsequence

  # Read, validate, and load in the FASTA data from path
  fastaDF <- loadFASTA(fastaPath)  # Dataframe for FASTA sequences

  # Initialize counters
  numMotif <- c(); lenSeq <- c()

  for (i in 1:length(fastaDF$seq)){

    # Iterate upstream and find any motifs that could be the
    #inputted binding motifs
    # (ie. find all matches of the inputted subsequence)
    motifIndex <- findMotifs(fastaDF$seq[i], pattern)
    mi <- attr(motifIndex, "match.length")

    # Below show zero when no match (implied by -1 match.length value)
    # else show number of matches
    numMotif[i] <- length(mi[which(mi != -1)])
    lenSeq[i] <- nchar(fastaDF$seq[i])
  }

  # How many motifs could fit each sequence at max
  lenpatt <- nchar(pattern) # Length (in nucleotides) of ortholog pattern

  if (lenpatt == 0){
    # No motif can be found if there is no motif to look for
    possibleMotif <- rep(0, length(fastaDF$seq))
  } else {
    # If entered motif is valid
    possibleMotif <- lenSeq/lenpatt # (Len of seq) / (Len of pattern)
  }

  # Create tibble values
  id <- seq(1, length(fastaDF$head))
  ratioMotif <- numMotif / possibleMotif

  # Create tibble
  motifDF <- tibble::tibble(id = id,
                            head = fastaDF$head,
                            numMotif = numMotif,
                            lenSeq = lenSeq,
                            ratioMotif = ratioMotif)

  return(motifDF) # Return the plot so it'll be savable outside function
}

#' Plot Ortholog Subsequence (gene/motif) Occurrence
#'
#' A function that plots how many times user-inputted ortholog pattern appears
#'    in each sequence in FASTA file. Each bar represents one sequence.
#'
#' @param fastaPath is a character string indicating the file path for input
#' FASTA file containing DNA sequences.
#' @param pattern is a character string indicating the ortholog gene/motif
#' subsequence to search for
#'
#' @return Returns a barplot of frequency (number of occurrence) for all
#'    sequences in FASTA file. Each bar represents a sequence, and bar height is
#'    the frequency (number of times pattern occurred) in that sequence.
#'
#' @examples
#'
#' #Example 1:
#' #Generate plot for FASTA file containing 1 DNA sequence (1 bar) and target
#' #orthologous binding motif pattern "CGCG".
#' # Get raw data
#' placenta_path <-  system.file("extdata", "P_Placenta_seq.fasta", package = "orthoVisualizer")
#' freqSeq(fastaPath = placenta_path, pattern = "CGCG")
#'
#' #Example 2:
#' #Generate plot for FASTA file containing 6 DNA sequences (6 bars) and target
#' #orthologous binding motif pattern "CGCG".
#' # Get raw data
#' DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
#' freqSeq(fastaPath = DNA_seq_path, pattern = "CGCG")
#'
#' @references
#' Bache, SM., et al. A Forward-Pipe Operator for R
#' R Package Magrittr Version 2.0.3. Comprehensive R Archive Network (CRAN),
#' 30 Mar. 2022, \href{https://cran.r-project.org/web/packages/magrittr/index.html}{Link}.
#'
#' Wickham, H., et al. A Grammar of Data Manipulation
#' R Package Dplyr Version 1.0.10. The Comprehensive R Archive Network (CRAN),
#' 1 Sept. 2022, \href{https://cran.r-project.org/web/packages/dplyr/index.html}{Link}.
#'
#' Wickham, H. et al. Create Elegant Data Visualisations Using the Grammar of
#' Graphics R Package GGPLOT2 Version 3.4.0.
#' The Comprehensive R Archive Network (CRAN), 4 Nov. 2022,
#' \href{https://cran.r-project.org/web/packages/ggplot2/index.html}{Link}.
#'
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
freqSeq <- function(fastaPath, pattern) {

  # Get quantity info
  motifDF <- quantSeq(fastaPath = fastaPath, pattern = pattern)

  # Title of plot
  plotTitle <- paste("Frequency (Number of Occurrence) of pattern",
                     pattern, "per Sequence")

  id <- numMotif <- NULL

  plot <- motifDF %>%
    ggplot2::ggplot(ggplot2::aes(x=id, y=numMotif))+
    ggplot2::geom_col(position = "dodge", fill="dark blue") +
    ggplot2::labs(title=plotTitle) +
    ggplot2::ylab("Frequency (Number of Occurence)") +
    ggplot2::xlab("i-th Sequence in input FASTA")

  plot # Display the plot

  return(plot) # Return the plot so it'll be savable outside function
}


#' Plot Ortholog Subsequence (gene/motif) Occurrence Ratio
#'
#' A function that plots ratio of user-inputted ortholog pattern appearance
#'    in each sequence in FASTA file. Each bar represents one sequence.
#'
#' @param fastaPath is a character string indicating the file path for input
#' FASTA file containing DNA sequences.
#' @param pattern is a character string indicating the ortholog gene/motif
#' subsequence to search for
#'
#' @return Returns a barplot of frequency ratio for all sequences in FASTA file.
#'    Each bar represents a sequence, and bar height is the frequency ratio
#'    in that sequence.
#'    Frequency ratio = (Number of Occurrence) / (Length of Seq / Length of pattern)
#'
#' @examples
#'
#' #Example 1:
#' #Generate plot for FASTA file containing 1 DNA sequence (1 bar) and target
#' #orthologous binding motif pattern "CGCG".
#' # Get raw data
#' placenta_path <-  system.file("extdata", "P_Placenta_seq.fasta", package = "orthoVisualizer")
#' freqRatioSeq(fastaPath = placenta_path, pattern = "CGCG")
#'
#' #Example 2:
#' #Generate plot for FASTA file containing 6 DNA sequences (6 bars) and target
#' #orthologous binding motif pattern "CGCG".
#' # Get raw data
#' DNA_seq_path <-  system.file("extdata", "DNA_seq.fasta", package = "orthoVisualizer")
#' freqRatioSeq(fastaPath = DNA_seq_path, pattern = "CGCG")
#'
#' @references
#' Bache, SM., et al. A Forward-Pipe Operator for R
#' R Package Magrittr Version 2.0.3. Comprehensive R Archive Network (CRAN),
#' 30 Mar. 2022, \href{https://cran.r-project.org/web/packages/magrittr/index.html}{Link}.
#'
#' Wickham, H., et al. A Grammar of Data Manipulation
#' R Package Dplyr Version 1.0.10. The Comprehensive R Archive Network (CRAN),
#' 1 Sept. 2022, \href{https://cran.r-project.org/web/packages/dplyr/index.html}{Link}.
#'
#' Wickham, H. et al. Create Elegant Data Visualisations Using the Grammar of
#' Graphics R Package GGPLOT2 Version 3.4.0.
#' The Comprehensive R Archive Network (CRAN), 4 Nov. 2022,
#' \href{https://cran.r-project.org/web/packages/ggplot2/index.html}{Link}.
#'
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
freqRatioSeq <- function(fastaPath, pattern) {

  # Get tibble containing quantified search results
  motifDF <- quantSeq(fastaPath = fastaPath, pattern = pattern)

  # Title of plot
  plotTitle <- paste("Frequency Ratio of pattern", pattern, "per Sequence")
  plotPosFreqCap <- paste("Possible Frequence = (Sequence length) / ",
                          "(pattern length = ",
                          nchar(pattern), ")", sep="")
  plotCaption <- paste("[Caption]: Frequency Ratio =",
                       "(Frequence) / (Possible Frequence)\n",
                       "Frequence = Actual number of occurrence of",
                       pattern,
                       "per sequence \n", plotPosFreqCap
  )

  id <- ratioMotif <- NULL

  plot <- motifDF %>%
    ggplot2::ggplot(ggplot2::aes(x=id, y=ratioMotif))+
    ggplot2::geom_col(position = "dodge", fill="dark red") +
    ggplot2::labs(title=plotTitle, caption = plotCaption) +
    ggplot2::ylab("Frequency Ratio") +
    ggplot2::xlab("i-th Sequence in input FASTA")

  plot # Display the plot

  return(plot) # Return the plot so it'll be savable outside function
}

# Create helper functions=====================

#' Loads FASTA files (helper function)
#' @param fastaPath Pathway of FASTA file with DNA Sequence
#'
#' @return Returns loaded FASTA in optimal format for analysis
#'
#' @references
#' Muller, K., et al. In-Line Documentation for R
#' R Package Roxygen2 Version 7.2.2. Comprehensive R Archive Network (CRAN),
#' 22 July. 2022, \href{https://cran.r-project.org/web/packages/roxygen2/index.html}{Link}.
#'
#' @importFrom tibble tibble
loadFASTA <- function(fastaPath){

  # Read FASTA file from its path
  fastaFile <- readLines(fastaPath)

  # Get location header appears
  lineheads <- grep("^[>]", fastaFile)

  # Insert headers
  fastahead <- fastaFile[lineheads]

  # Validate sequences
  if (length(grep("[^ATCG]", fastaFile[-lineheads]))!= 0) {
    stop("ERROR: FASTA formatted incorrectly: nucleotide sequence contains newline or non-ATCG characters.")
  }

  ithSeq <- c(); start <- c(); end <- c();

  # Insert corresponding sequences
  for (i in 1:length(lineheads)){
    if (i == length(lineheads)){
      start[i] <- (lineheads[i] + 1)
      end[i] <- length(fastaFile)
      ithSeq[i] <- paste(fastaFile[start[i]:end[i]], collapse = "", sep="")
    } else {
      start[i] <- (lineheads[i] + 1)
      end[i] <- (lineheads[i+1] - 1)
      ithSeq[i] <- paste(fastaFile[start[i]:end[i]], collapse = "", sep="")
    }
  }  # End of for loop of i

  # Create tibble
  fastaDF <- tibble::tibble(head = fastahead, seq = ithSeq, start=start, end=end)

  # Return the dataframe
  return(fastaDF)
}


findMotifs <- function(seq, pattern) {

  patt <- paste(".",pattern,".", sep="")  # Get regex expression ".[pattern]."
  matches <- gregexpr(patt, seq)  # Return start index of each match found

  # This vector contains starting indices of all matches w/o names
  # Simply a vector containing indicies but we do not know which substring
  # this starting index corresponds to
  motifsIndex <- matches[[1]]

  # Store all matched substrings and treat them as names of motifIndex vector
  # i.e. Connect matched motifs (substrings) to its corresponding starting index
  names(motifsIndex) <- regmatches(seq, matches)[[1]]

  # Contains indices and actual substrings that matched Regex
  return(motifsIndex)
}

#' Print Annotation of Ortholog Subsequence (gene/motif) Occurrence
#'
#' @param header is the header of sequence
#' @param seq is the DNA sequence to print annotation
#' @param pattern is a character string indicating the ortholog gene/motif
#' subsequence to search for .
#' @param motifIndex is list of indices of where the pattern is found in seq
#'
#' @return Returns invisible NULL. Prints each sequence in the FASTA file with
#'    annotation.
#'
printAnnotate <- function(header, seq, pattern, motifIndex){

  # Print header at the very top
  cat(header, '\n')

  # Print 5'- : Start of sequence
  cat("5'-")

  # Set indent for aligning the sequences
  indent <- "   "  # Indent for formatting

  # Print seq genome sequence w/ underline annotation
  for (i in 1:length(seq)) {

    # Print lines of nucleotides
    if (i != 1) {
      # If not first line, print indent to align seq.
      cat(indent)
    }

    # Print i-th line of 50 nt in sequence
    if (i == length(seq)){
      cat(seq[i], "-3'\n", sep="")
    } else{
      cat(seq[i], '\n')
    }

    # Create string containing ALL annotations for i-th line;
    # Currently 50 whitespaces, where some of them will be replaced with motif
    # annotations.
    annotations <- strrep(" ", 50)

    # Find motif substrings with its start index on line i
    sel <- ((50*i - 49) <= motifIndex) & (motifIndex <= 50*i)  # Find w/ indices
    # Find names (i.e Actual 6-char substrings of pattern ".CGCG.")
    motifs <- names(motifIndex[sel])  # Vector with all motifs in line i

    patt_sign <- rep("=", nchar(pattern))
    patt_sign <- paste(patt_sign, sep="", collapse = "")

    # Change the motif names to motif annotations by substitution
    motifs_ann <- gsub(pattern, patt_sign, motifs)# Change middle 'GCGC' to '==='

    # Know the motif length so that we know where is the end index (ie. end = start + (motif length) - 1)
    motifLength <- attr(motifIndex, "match.length")[1]

    # Combine all motif annotations for line i into single string
    # i.e) Select the blanks and replace corresponding indices with the
    # annotation.
    for (j in seq_along(motifs_ann)) {
      mI <- motifIndex[sel][j]   # Starting index of the motif we are looking at
      start <- mI - 50*(i-1)     # Position (1-50) the motif starts on line i
      end <- start + motifLength - 1   # Index where motifLength-char-long motif ends on line i

      # The white spaces between the start~end index is replaced with the
      # corresponding motifLength-char motif annotation.
      substring(annotations, start, end) <- motifs_ann[j]
    }

    # Print annotations for line i
    cat(indent, annotations, '\n', sep="")


  } # End of for-loop over seq.

  return(invisible(NULL))
}
