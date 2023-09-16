# Metagenomics
This application takes a file containing reads and maps them to the genomes provided. The complexity/difficulty lies in that the reads can originate from any of the 1000 genomes provided. As a result, one must implement read mapping capabilities for one genome and then apply that solution across all of the genomes efficiently.

## Input Files
- reads.fasta - file containing 20,000 reads, each with a unique ID
- genome\_{ID}.fasta - 1,000 reference genomes

## Output Files
- predictions.csv - resulting file that maps each read to the genome it originated from. The format is '>read\_{READ_ID} Genome_Number{GENOME_ID}'. Note, the desired output is to map each read to a genome so all reads must be mapped to a genome (predictions.csv should be 20,000 reads).

## Usage
After setting the desired constants for the specified use case and data, this can be run with `python3 main.py`. Sample reads and sample genomesare provided in the data directory.

## Results
- Successfully identifies present genomes based on custom threshold.
- Performs metagenomics analysis to map each read to a genome, account for SNPs in the form of substitutions.
- Successfully mapped 90% of reads to the correct genome in case of 20,000 reads with 1,000 possible genomes each of length 10,000 in under a minute.

### Todo / Future Improvements
- Apply further optimizations such as bloom filters, minimizers, or minimum perfect hashing to enable scaling for even larger datasets.
- Increase accuracy in detecting indels and duplicated regions between genomes with respect to mapping a read to a genome.
