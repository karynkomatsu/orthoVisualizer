
# OrthoVisualizer -- package for genome analysis (orthologue identification and visualization) of various species


# Create helper functions=====================

loadFASTA <- function(fastaPath){
  # Input: pathway of FASTA file

  # Read and validate FASTA. Then create and return a tibble with columns:
  #    Header: Header of the FASTA file
  #    Sequence: Entire sequence in the FASTA
  #    Upstream Sequence: First (n-k) nucleotides of the FASTA file
  #    Coding Sequence: Last k/3 codons (k nucleotides) of the FASTA file

  # Read FASTA file from its path
  fastaFile <- readLines(fastaPath)

  # TODO Delete this line later
  #fastaFile <- readLines("Data/assembly.fasta")#[1:8606]

  # Get location header appears
  lineheads <- grep("^[>]", fastaFile)

  # Insert headers
  fastahead <- fastaFile[lineheads]

  # Validate sequences
  if (length(grep("[^ATCG]", fastaFile[-lineheads]))!= 0) {
    print("ERROR: FASTA formatted incorrectly: nucleotide sequence contains newline or non-ATCG characters.")
  }

  ithSeq <- c(); start <- c(); end <- c();

  # Insert corresponding sequences
  for (i in 1:length(lineheads)){
    if (i == length(lineheads)){
      start[i] <- (lineheads[i] + 1)
      end[i] <- length(fastaFile)
      ithSeq[i] <- BiocGenerics::paste(fastaFile[start[i]:end[i]], collapse = "", sep="")
    } else {
      start[i] <- (lineheads[i] + 1)
      end[i] <- (lineheads[i+1] - 1)
      ithSeq[i] <- BiocGenerics::paste(fastaFile[start[i]:end[i]], collapse = "", sep="")
    }
  }  # End of for loop of i

  # Create tibble
  # head = header; seq = corresponding sequence; start = line that seq started; end = line that seq ended
  fastaDF <- tibble::tibble(head = fastahead, seq = ithSeq, start=start, end=end)

  # Separate the sequence into "upstream" and "coding" sequences
  # ie. We add 2 more columns to tibble of fasta returned by readFASTA()
  fastaDF$upstream <- fastaDF$seq #substring(fastaDF$seq, 1, 500)
  fastaDF$coding <- substring(fastaDF$seq, length(fastaDF$seq) - 30, length(fastaDF$seq)) #substring(fastaDF$seq, 501, 530)

  # Check if sequence splitted
  #if (paste(fastaDF$upstream, fastaDF$coding, sep="")!= paste(fastaDF$seq)) {
  #  print("ERROR: FAILED to split the sequence into upstream and coding...")
  #}

  # Return the edited dataframe with 2 extra columns
  return(fastaDF)
}


findMotifs <- function(upstream, patt = "CGCG") {
  # Input string: upstream of coding sequence

  # Iterate upstream sequence to find and return vector of indices that might be
  # motif binding regions (ex. "CGCG" pattern by default).

  patt <- paste(".",patt,".", sep="")  # Get regex expression ".[pattern]."
  matches <- gregexpr(patt, upstream)  # Return start index of each match found

  # This vector contains starting indices of all matches but DOES NOT HAVE NAMES
  # It is simply a vector containing indicies but we do not know which substring
  # this starting index corresponds to
  motifsIndex <- matches[[1]]

  # Store all matched substrings and treat them as names of motifIndex vector
  # i.e. Connect matched motifs (substrings) to its corresponding starting index
  names(motifsIndex) <- regmatches(upstream, matches)[[1]]

  # Contains indices and actual substrings that matched Regex
  return(motifsIndex)
}


getAA <- function(coding){
  # Input string: Coding sequence (should be 10 codon for this Unit)
  # Translate codons and return amino acid sequence with corresponding codons.

  # Validate input: Is the input string divisible by 3? (ie NO incomplete codon)
  if (nchar(coding) %% 3 != 0){
    print("INVALID INPUT: Coding sequence string is non-splittable by 3.")
    return(invisible(NULL))
  }

  # GENETIC CODE: Linking Codons to Amino Acids
  code <- Biostrings::GENETIC_CODE

  # Translate individual codons to AA using loop
  start <- seq(1, nchar(coding), 3)   # starting indices of codons: 1, 4, 7,...
  # Vector that will contain AA
  AAseq <- c()

  # Loop starts
  for (i in 1:length(start)){
    codon <- substring(coding, start[i], start[i] + 2)
    AAseq[i] <- code[codon]    # Add translated AA
    names(AAseq)[i] <- codon   # Add corresponding codon as name
  }

  return(AAseq)
}


printAnnotate <- function(header, upstream, pattern, motifIndex){
  # Inputs:
  #    header: Header of the FASTA file
  #    upstream: Vector of strings representing upstream region.
  #              Each item in vector is 50 nucleotide long string.
  #    motifIndex: Vector of motifs' starting indices associated with their
  #               motif substrings.
  #    AAseq: Sequence of translated amino acids associated with their codons

  # Format our annotation results as shown in the "Sample Annotation" section
  # under Professor Steipe's Integrator Unit: Genome Annotation.

  # Print header at the very top
  cat(header, '\n')

  # Print 5'- : Start of sequence
  cat("5'-")

  # Set indent for aligning the sequences
  indent <- "   "  # We add this indent to format the upstream

  # Print upstream genome sequence w/ underline annotation
  for (i in 1:length(upstream)) {

    # Print lines of nucleotides
    if (i != 1) {
      # If not first line, print indent to align seq.
      cat(indent)
    }

    # Print i-th line of 50 nt in upstream sequence
    if (i == length(upstream)){
      cat(upstream[i], "-3'\n", sep="")
    } else{
      cat(upstream[i], '\n')
    }

    # Create string containing ALL annotations for i-th line;
    # Currently 50 whitespaces, where some of them will be replaced with motif
    # annotations.
    annotations <- strrep(" ", 50)

    # Find motif substrings with its start index on line i
    sel <- ((50*i - 49) <= motifIndex) & (motifIndex <= 50*i)  # Find w/ indices
    # Find names (i.e Actual 6-char substrings of pattern ".CGCG.")
    motifs <- names(motifIndex[sel])  # Vector with all motifs in line i

    # Change the motif names to motif annotations by substitution
    ## TODO: UNcomment the motifs_ann below
    motifs_ann <- gsub(pattern, "====", motifs)# Change middle 'GCGC' to '==='
    #motifs_ann <- gsub("[AT]", "=", motifs_ann)   # change A or T at end to '='
    #motifs_ann <- gsub("[CG]", "X", motifs_ann)   # Change C or G at end to 'X'

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


  } # End of for-loop over upstream seq.

  return(invisible(NULL))
}


# Function that combines all the helper functions================

# For sequence annotation of specific subsequence pattern representing otholog gene/motif
annotateSeq <- function(fastaPath="./inst/extdata/assembly.fasta", pattern="CGCG") {

  # Input: Path to FASTA file, Pattern of the target orthologous motif / gene subsequence

  # Annotate input FASTA sequence with Mbp1 binding motifs potentially observed.
  # Gathering required information from the FASTA file and combine the helper
  # functions to print the annotated result.
  # Print: Formatted annotation of FASTA file at the inputted file path.

  # Check if required packages installed
  if (! requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (! requireNamespace("Biostrings", quietly = TRUE)) {
    BiocManager::install("Biostrings")
  }

  # Read, validate, and load in the FASTA data from path
  fastaDF <- loadFASTA(fastaPath)  # Dataframe for FASTA sequences

  for (i in 1:length(fastaDF$seq)){
    # Get which of  split the upstream into lines of 50 nucleotides
    start <- seq(1, nchar(fastaDF$seq[i]), 50)             # start indices
    stop <- seq(50, nchar(fastaDF$seq[i]), 50)             # stop indices

    # Add the final stop index AKA the final incomplete line
    if (!(nchar(fastaDF$seq[i]) %in% stop)){
      stop <- append(stop, nchar(fastaDF$seq[i]))
    }

    # Add all the lines of sequences
    splitUpstream <- substring(fastaDF$upstream, start, stop)

    # Iterate upstream and find any motifs that could be the inputted binding motifs
    # (ie. find all matches of the inputted subsequence)
    # Vectors of motifs & starting indices
    motifIndex <- findMotifs(fastaDF$seq[i], pattern)

    # Print formatted genome sequence with corresponding annotations
    printAnnotate(fastaDF$head[i], splitUpstream, pattern, motifIndex)
  }
  return(invisible(NULL))
}


# frequency of a sequence
freqSeq <- function(fastaPath="./inst/extdata/assembly.fasta", pattern="CGCG") {

  # Input: Path to FASTA file,
  #        Pattern of the target orthologous motif / gene subsequence

  # Check if required packages installed
  if (! requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (! requireNamespace("Biostrings", quietly = TRUE)) {
    BiocManager::install("Biostrings")
  }
  if (! requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (! requireNamespace("magrittr", quietly = TRUE)) {
    install.packages("magrittr")
  }
  if (! requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (! requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse")
  }

  loadNamespace("ggplot2")
  loadNamespace("tidyverse")
  loadNamespace("magrittr")
  loadNamespace("dplyr")

  # Read, validate, and load in the FASTA data from path
  fastaDF <- loadFASTA(fastaPath)  # Dataframe for FASTA sequences

  # Initialize counter for number of matching motifs per sequence
  numMotif <- c()

  for (i in 1:length(fastaDF$seq)){

    # Iterate upstream and find any motifs that could be the
    #inputted binding motifs
    # (ie. find all matches of the inputted subsequence)
    # Vectors of motifs & starting indices
    motifIndex <- findMotifs(fastaDF$seq[i], pattern)
    mi <- attr(motifIndex, "match.length")
    numMotif[i] <- length(mi[which(mi != -1)])

  }
  # Create tibble
  motifDF <- tibble::tibble(id = seq(1, length(fastaDF$head)),
                            head = fastaDF$head,
                            numMotif = numMotif)

  # Title of plot
  plotTitle <- paste("Frequency (Number of Occurrence) of pattern",
                     pattern, "per Sequence")

  plot <- motifDF %>%
    ggplot2::ggplot(ggplot2::aes(x=id, y=numMotif))+
    ggplot2::geom_col(position = "dodge", fill="dark blue") +
    ggplot2::labs(title=plotTitle) +
    ggplot2::ylab("Frequency (Number of Occurence)") +
    ggplot2::xlab("i-th Sequence in input FASTA")

  plot # Display the plot

  return(plot) # Return the plot so it'll be savable outside function
}

# frequency RATIO of a sequence
freqRatioSeq <- function(fastaPath="./inst/extdata/assembly.fasta", pattern="CGCG") {

  # Input: Path to FASTA file,
  #        Pattern of the target orthologous motif / gene subsequence

  # Check if required packages installed
  if (! requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (! requireNamespace("magrittr", quietly = TRUE)) {
    install.packages("magrittr")
  }
  if (! requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (! requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse")
  }

  loadNamespace("ggplot2")
  loadNamespace("tidyverse")
  loadNamespace("magrittr")
  loadNamespace("dplyr")

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
                            numMotif = numMotif,
                            lenSeq = lenSeq,
                            ratioMotif = ratioMotif)

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

  plot <- motifDF %>%
    ggplot2::ggplot(ggplot2::aes(x=id, y=ratioMotif))+
    ggplot2::geom_col(position = "dodge", fill="dark red") +
    ggplot2::labs(title=plotTitle, caption = plotCaption) +
    ggplot2::ylab("Frequency Ratio") +
    ggplot2::xlab("i-th Sequence in input FASTA")

  plot # Display the plot

  return(plot) # Return the plot so it'll be savable outside function
}

# Calling the Functions Above ====================

#annotateSeq(fastaPath="./inst/extdata/assembly.fasta")
#annotateSeq(fastaPath="./inst/extdata/POSPL_seq.fasta")
