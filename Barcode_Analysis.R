# BIOL432 Assignment #4, written by Zoe Clarke, submitted Jan 31, 2019

library(sangerseqR)
library(Biostrings)

# Read the BarcodePlateStats data file
stats <- read.csv("./Data/BarcodePlateStats.csv")

# Make an empty character vector with the length equal to the number of
# Ok DNA sequences
goodFiles <- vector(mode = "character", length = sum(stats$Ok == TRUE))

# Make a loop that goes through every row in "stats" and puts the names of the sequences
# that are Ok in a list
j = 1
for (i in 1:length(stats$Ok)) {
  if (stats$Ok[i] == TRUE) {
    goodFiles[j] <- as.character(stats$Chromatogram[i])
    j = j + 1
  }
}

# Create an empty vector to put the fasta formatted sequences in
sequences <- vector(mode = "character")

# Create an empty vector to put the file names into
header <- vector(mode = "character")

# This number will be used to assign the sequence to the vector
i = 1

# Below is the loop to read the files and put the primary sequences into a vector

for (file in goodFiles) {
  ITS <- read.abif(paste("./Data/", file, sep = "")) # Read file from Data folder
  ITSseq <- sangerseq(ITS) # Extract file
  SeqX <- makeBaseCalls(ITSseq) # Call file
  primary <- primarySeq(SeqX, string = TRUE) # Get primary sequence only in character type
  header[i] <- file # Put file name into a vector
  sequences[i] <- primary # Put sequence into a vector
  i = i + 1
}

# Load a package to make fasta files; warning, do not load earlier!
library(seqinr)

# Use a loop to place each sequence into the same .fasta file called "output"
for (y in 1:length(header)) {
  # Makes the output file with the first sequence
  if (y == 1) {
    write.fasta(sequences[y], file.out = "output", as.string = TRUE,
                names = paste(header[y], " length=", nchar(sequences[y]), 
                              "; type=dna", sep = ""),
                open = "w")
  }
  # Appends the rest of the sequences and headers onto the output.fasta file
  else if (y > 1) {
    write.fasta(sequences[y], file.out = "output", as.string = TRUE,
                names = paste(header[y], " length=", nchar(sequences[y]), 
                              "; type=dna", sep = ""),
                open = "a")
      }
}

# You should have an output.fasta file of the primary sequences