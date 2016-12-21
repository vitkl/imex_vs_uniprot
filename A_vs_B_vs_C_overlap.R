## Function which compares A with B from logic table
## Generates summary


A_vs_B_vs_C_overlap =  function(proteome_vs_interactome_o, A, B, C, SPECIES_NAME, reviewed, isoforms) {

  A_vs_B_vs_C = data.frame(
  "species_name" = SPECIES_NAME,
  "A" = A,
  "B" = B,
  "C" = C,
  "reviewed" = reviewed,
  "isoforms" = isoforms,
  "count_A" = sum((proteome_vs_interactome_o[A] == 1)),
  "count_B" = sum((proteome_vs_interactome_o[B] == 1)),
  "count_C" = sum((proteome_vs_interactome_o[C] == 1)),
  "A_exclusively" = sum((proteome_vs_interactome_o[A] == 1) & (0 == proteome_vs_interactome_o[B]) & (0 == proteome_vs_interactome_o[C])),
  "A_and_B"       = sum((proteome_vs_interactome_o[A] == 1) & (1 == proteome_vs_interactome_o[B])),
  "A_and_B_not_C" = sum((proteome_vs_interactome_o[A] == 1) & (1 == proteome_vs_interactome_o[B]) & (0 == proteome_vs_interactome_o[C])),
  "B_exclusively" = sum((proteome_vs_interactome_o[A] == 0) & (1 == proteome_vs_interactome_o[B]) & (0 == proteome_vs_interactome_o[C])),
  "A_and_C"       = sum((proteome_vs_interactome_o[A] == 1) & (1 == proteome_vs_interactome_o[C])),
  "A_and_C_not_B" = sum((proteome_vs_interactome_o[A] == 1) & (0 == proteome_vs_interactome_o[B]) & (1 == proteome_vs_interactome_o[C])),
  "C_exclusively" = sum((proteome_vs_interactome_o[A] == 0) & (0 == proteome_vs_interactome_o[B]) & (1 == proteome_vs_interactome_o[C])),
  "B_and_C"       = sum((proteome_vs_interactome_o[B] == 1) & (1 == proteome_vs_interactome_o[C])),
  "B_and_C_not_A" = sum((proteome_vs_interactome_o[A] == 0) & (1 == proteome_vs_interactome_o[B]) & (1 == proteome_vs_interactome_o[C])),
  "A_and_B_and_C" = sum((proteome_vs_interactome_o[A] == 1) & (1 == proteome_vs_interactome_o[B]) & (1 == proteome_vs_interactome_o[C])),
  stringsAsFactors = FALSE
)
colnames(A_vs_B_vs_C) = c("species_name", 
                     "A", 
                     "B", 
                     "C",
                     "reviewed", 
                     "isoforms", 
                     paste("count_", A, sep = ""), 
                     paste("count_", B, sep = ""), 
                     paste("count_", C, sep = ""),
                     paste(A, "_exclusively", sep = ""),
                     paste(A, "_and_", B, sep = ""),
                     paste(A, "_and_", B, "_not_", C, sep = ""), 
                     paste(B, "_exclusively", sep = ""),
                     paste(A, "_and_", C, sep = ""),
                     paste(A, "_and_", C, "_not_", B, sep = ""),
                     paste(C, "_exclusively", sep = ""),
                     paste(B, "_and_", C, sep = ""),
                     paste(B, "_and_", C, "_not_", A, sep = ""),
                     paste(A, "_and_", B, "_and_", C, sep = "")
                     )

return(A_vs_B_vs_C)

}