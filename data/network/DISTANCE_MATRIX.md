# Enzyme-metabolite distance matrix

`model_dist.bin` (the full Human1 distance matrix) is **720,784,512 bytes** =
9492 x 9492 IEEE-754 little-endian doubles, **column-major** (the order produced by
R's `as.vector()`), and is therefore distributed outside Git:

    SHA-256  bf8d681b2d0c8065f5bde01ce96d656eca0ca243ed78fdf0fcac1a0c2a65a093

Rebuild it from the source model with

    Rscript scripts/export_model.R        # writes model_dist.bin, enzyme.csv, nodes.csv

`scripts/export_model.R` reads BioModels **MODEL2204280001 (Human1)**, converts the
reaction triples to an igraph object, and writes the shortest-path distance matrix in
which directly connected metabolite-enzyme pairs carry distance 2 (twice the hop count).
The three table files that ship with this repository (`nodes.csv`, `enzyme.csv`,
`nodeid_map.csv`, `EC_gene_map.csv`) are the inputs and outputs of that step and are
checksummed in `CHECKSUMS.sha256`.

Verify a rebuilt matrix against the archived one before using it:

    sha256sum model_dist.bin
