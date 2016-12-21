## Function which compares A with B from logic table
## Generates summary


A_vs_B_overlap =  function(proteome_vs_interactome_o, A, B, SPECIES_NAME, reviewed, isoforms) {

A_vs_B = data.frame(
  "species_name" = SPECIES_NAME,
  "A" = A,
  "B" = B,
  "reviewed" = reviewed,
  "isoforms" = isoforms,
  "count_A" = sum((proteome_vs_interactome_o[A] == 1)),
  "count_B" = sum((proteome_vs_interactome_o[B] == 1)),
  "A_exclusively" = sum((proteome_vs_interactome_o[A] == 1) & (0 == proteome_vs_interactome_o[B])),
  "A_and_B" = sum((proteome_vs_interactome_o[A] == 1) & (1 == proteome_vs_interactome_o[B])),
  "B_exclusively" = sum((proteome_vs_interactome_o[A] == 0) & (1 == proteome_vs_interactome_o[B]))
)
colnames(A_vs_B) = c("species_name", "A", "B", "reviewed", "isoforms", 
                     paste("count_", A, sep = ""), paste("count_", B, sep = ""), 
                     paste(A, "_exclusively", sep = ""),
                     paste(A, "_and_", B, sep = ""), paste(B, "_exclusively", sep = ""))

return(A_vs_B)

}