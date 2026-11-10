library(digest)

# ANI from Jaccard
compute_ANI <- function(k, J) {
  if (J <= 0) return(NA)
  ANI <- 1 + (1/k) * log((2 * J) / (1 + J))
  return(ANI)
}

# Generate random sequence
generate_sequence <- function(length) {
  sample(c("A", "C", "G", "T"), length, replace = TRUE)
}

# Extract k-mers from sequence
extract_kmers <- function(seq, k) {
  if (length(seq) < k) return(character(0))
  kmers <- character(length(seq) - k + 1)
  for (i in 1:(length(seq) - k + 1)) {
    kmers[i] <- paste0(seq[i:(i+k-1)], collapse = "")
  }
  return(kmers)
}

# Hash k-mers to integers for minHash
hash_kmers <- function(kmers) {
  if (length(kmers) == 0) return(numeric(0))
  sapply(kmers, digest, algo = "md5", serialize = FALSE)
}

# MinHash sketch: keep k smallest hashes
sketch <- function(kmers, sketch.size) {
  if (length(kmers) == 0) return(numeric(0))
  hashes <- hash_kmers(kmers)
  sketch <- sample(hashes, sketch.size)
  return(sketch)
}

# Calculate Jaccard from two sketches
jaccard_from_sketches <- function(sketch1, sketch2) {
  if (length(sketch1) == 0 || length(sketch2) == 0) return(0)
  intersection <- length(intersect(sketch1, sketch2))
  union <- length(union(sketch1, sketch2))
  if (union == 0) return(0)
  return(intersection / union)
}

# k-mer size
k <- 21

# Generate reference genome
genome.size <- 5e6
genome.ref <- generate_sequence(genome.size)
kmers.ref <- extract_kmers(genome.ref, k)

sketch.list <- c(500, 1000, 2000, 5000, 10000)
props.list <- seq(from = 0.7, to = 1.0, by = 0.05)

# Create results list with one data frame per k-mer size
results_list <- list()

# Create matrix for this k value
results_matrix <- matrix(NA, 
                        nrow = length(sketch.list), 
                        ncol = length(props.list))
colnames(results_matrix) <- paste0("p", props.list)
rownames(results_matrix) <- paste0("sk", sketch.list)

# Simulation for this k value
p <- 1
for (p in 1:length(props.list)) {
  prop.match <- props.list[p]
  
  # Generate query genome with prop.match sites identical to reference
  genome.query <- genome.ref
  n_mismatches <- round((1 - prop.match) * genome.size)
  if (n_mismatches > 0) {
    mismatch_positions <- sample.int(genome.size, n_mismatches)
    # Ensure mutated base is different from original
    for (i in 1:n_mismatches) {
      pos <- mismatch_positions[i]
      orig_base <- genome.ref[pos]
      possible_bases <- c("A", "C", "G", "T")[c("A", "C", "G", "T") != orig_base]
      genome.query[pos] <- sample(possible_bases, 1)
    }
  }
  
  # Extract k-mers
  kmers.query <- extract_kmers(genome.query, k)
  
  # Calculate true Jaccard (no sketching)
  true_jaccard <- length(intersect(kmers.ref, kmers.query)) / length(union(kmers.ref, kmers.query))
  
  # Loop over sketch sizes
  s <- 1
  for (s in 1:length(sketch.list)) {
    sketch.size <- sketch.list[s]
    
    # Create MinHash sketches
    sketch.ref <- minhash_sketch(kmers.ref, sketch.size)
    sketch.query <- minhash_sketch(kmers.query, sketch.size)
    
    # Calculate Jaccard from sketches
    J <- jaccard_from_sketches(sketch.ref, sketch.query)
    
    # Store Jaccard (not ANI, to see raw relationship)
    results_matrix[s, p] <- compute_ANI(k, J)
  }
}

results_matrix
