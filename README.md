# CTCFResearch

**Immediate Goal:**

- Start with Chromosome 21 (chr21) because it is the smallest chromosome and will be easier to test and debug before scaling to the full dataset.

**Understand the BED File Format**
- Obtain and inspect the BED files.
- Each line appears to represent an edge/connection in a graph.
- Expected fields are roughly:
  - Chromosome/coordinate information for one genomic location
  - Chromosome/coordinate information for another genomic location
  - A value representing the strength/weight of the relationship
- Treat the final relationship-strength value as the edge weight in the graph.

**Design the Graph Representation**
- Determine how genomic regions will map to graph nodes.
- Determine how BED file rows will map to graph edges.
- Build an efficient graph structure:
  - Prefer an adjacency list for scalability.
  - Consider how to generate an adjacency matrix when needed.
- Focus on efficient storage and retrieval so processing larger chromosomes does not become prohibitively slow. 

**Build the Initial Pipeline**
1. Load BED file data. 
2. Parse connections for chr21 only.
3. Construct the graph. 
4. Verify that edge weights are correctly assigned. 
5. Test runtime and memory usage on the small chromosome before expanding. 
  
**Graph Analysis**
Calculate graph/network statistics such as:
- Density
- Degree distribution
- Connectivity metrics
- Other general graph statistics that may indicate meaningful structure or significance
  (Note: "statistics" here refers to graph/network statistics, not classical statistical tests or p-values.)
  
**Scale-Up Plan**
1. Get the BED files. 
2. Parse chr21. 
3. Build and validate graph representation. 
4. Compute graph statistics. 
5. Optimize performance if needed. 
6. Expand to additional chromosomes. 
7. Eventually process the full genome dataset.