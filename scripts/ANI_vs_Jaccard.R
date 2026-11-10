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

# Extract k-mers from sequence (vectorized)
extract_kmers <- function(seq, k) {
  if (length(seq) < k) return(character(0))
  n <- length(seq) - k + 1
  # Use sapply with window indices for faster k-mer extraction
  sapply(0:(n-1), function(i) {
    paste0(seq[(i+1):(i+k)], collapse = "")
  })
}

# Hash k-mers to integers for minHash (using MD5 converted to numeric)
hash_kmers <- function(kmers) {
  if (length(kmers) == 0) return(numeric(0))
  # Use MD5 hash converted to numeric (take first 8 hex chars)
  sapply(kmers, function(x) {
    h <- digest(x, algo = "md5", serialize = FALSE)
    as.numeric(paste0("0x", substr(h, 1, 8)))
  })
}

# MinHash sketch: keep k smallest hashes (vectorized)
minhash_sketch <- function(kmers, sketch.size) {
  if (length(kmers) == 0) return(numeric(0))
  hashes <- hash_kmers(kmers)
  if (length(hashes) <= sketch.size) return(sort(hashes))
  sort(hashes)[1:sketch.size]
}

# Calculate Jaccard from two sketches
jaccard_from_sketches <- function(sketch1, sketch2) {
  if (length(sketch1) == 0 || length(sketch2) == 0) return(0)
  intersection <- length(intersect(sketch1, sketch2))
  union <- length(union(sketch1, sketch2))
  if (union == 0) return(0)
  return(intersection / union)
}

# Mutate genome (vectorized)
mutate_genome <- function(genome.ref, prop.match) {
  genome.query <- genome.ref
  n_mismatches <- round((1 - prop.match) * length(genome.ref))
  
  if (n_mismatches > 0) {
    mismatch_positions <- sample.int(length(genome.ref), n_mismatches)
    # Vectorized mutation: sample from 3 possible bases
    bases <- c("A", "C", "G", "T")
    for (i in seq_along(mismatch_positions)) {
      pos <- mismatch_positions[i]
      genome.query[pos] <- sample(bases[bases != genome.ref[pos]], 1)
    }
  }
  return(genome.query)
}

# k-mer size
k <- 21

# Generate reference genome once
genome.size <- 5e7
genome.ref <- generate_sequence(genome.size)
kmers.ref <- extract_kmers(genome.ref, k)

sketch.list <- c(10000, 20000, 30000)
props.list <- seq(from = 0.7, to = 1.0, by = 0.05)

sketch.ref <- minhash_sketch(kmers.ref, max(sketch.list))  # Pre-compute sketch for reference

# Results matrix
results_matrix <- matrix(NA, 
                        nrow = length(sketch.list), 
                        ncol = length(props.list))
colnames(results_matrix) <- paste0("p", props.list)
rownames(results_matrix) <- paste0("sk", sketch.list)

# Simulation
for (p in seq_along(props.list)) {
  print(paste("P:", p))
  prop.match <- props.list[p]
  
  # Generate query genome with mutations
  genome.query <- mutate_genome(genome.ref, prop.match)
  kmers.query <- extract_kmers(genome.query, k)
  
  # Create query sketch
  sketch.query <- minhash_sketch(kmers.query, max(sketch.list))
  
  # Loop over sketch sizes
  for (s in seq_along(sketch.list)) {
    print(paste("S:", s))
    sketch.size <- sketch.list[s]
    
    # Use pre-computed reference sketch (subset to sketch.size)
    sketch.ref.subset <- sort(sketch.ref)[1:min(sketch.size, length(sketch.ref))]
    sketch.query.subset <- sort(sketch.query)[1:min(sketch.size, length(sketch.query))]
    
    # Calculate Jaccard from sketches
    J <- jaccard_from_sketches(sketch.ref.subset, sketch.query.subset)
    
    # Store ANI
    results_matrix[s, p] <- compute_ANI(k, J)
  }
}

results_matrix
