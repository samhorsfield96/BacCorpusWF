# ANI from Jaccard (corrected parentheses)
compute_ANI <- function(k, J) {
  ANI <- 1 + (1/k) * log((2 * J) / (1 + J))
  return(ANI)
}

genome.size <- 3e6
genome.ref <- rep(1, genome.size)

k <- 21

# sketch sizes to test
sketch.list <- c(1000, 2000, 3000, 4000, 5000,
                 6000, 7000, 8000, 9000, 10000)

# proportional overlap between genomes
props.list <- seq(from = 0.1, to = 1.0, by = 0.1)

# results matrix: rows = sketch sizes, cols = proportion match
results <- matrix(NA, 
                  nrow = length(sketch.list), 
                  ncol = length(props.list))

colnames(results) <- paste0("p", props.list)
rownames(results) <- paste0("sk", sketch.list)

# simulation
p <- 1
for (p in 1:length(props.list)) {
  
  prop.match <- props.list[p]
  
  n_match_positions <- round(prop.match * genome.size)
  
  # choose matching positions (these positions have identical k-mers in ref & query)
  match_positions <- sample.int(genome.size, n_match_positions)
  
  # create query genome with given proportion of shared k-mers
  genome.ref <- integer(genome.size)
  genome.query <- integer(genome.size)
  genome.ref[match_positions] <- 1
  genome.query[match_positions] <- 1
  
  # loop over sketch sizes
  s <- 1
  for (s in 1:length(sketch.list)) {
    
    sketch.size <- sketch.list[s]
    
    sk_ref_pos   <- sample.int(genome.size, sketch.size, replace = FALSE)
    sk_query_pos <- sample.int(genome.size, sketch.size, replace = FALSE)
    
    common_positions <- intersect(sk_ref_pos, sk_query_pos)
    intersection_count <- sum(common_positions %in% match_positions)
    
    # sample sketches

    union_positions <- union(sk_ref_pos[sk_ref_pos %in% match_positions],
                             sk_query_pos[sk_query_pos %in% match_positions])
    union_count <- length(union_positions)
    
    # handle degenerate case: if union_count == 0, set J = 0 (no matching k-mers sampled)
    if (union_count == 0) {
      J <- 0
    } else {
      J <- intersection_count / union_count
    }
    
    # compute ANI from Jaccard
    ANI.est <- compute_ANI(k, J)
    
    results[s, p] <- ANI.est
  }
}

results
