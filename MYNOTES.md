**Implementation Reasoning**

My first concern regarding implementation was about the magnitude of the files. 
These .bed files containing chromosome content could be 1000s of lines long, and the 
usual dataframe or Networkxx could not be used efficiently. I decided 
to attack this project with a new implementation: Parquets.

'.parquet' allows data to be stored in columns as compact 
binary numbers, with built-in statistics and compression. Plain text (.tsv, .csv, etc) 
have every number stored as human-readable, which can slow the runtime significantly. 
Parquet can skip over irrelevant chunks/columns instead of parsing every single line 
from scratch.

**Tools Used:**

1. `polars` + `pyarrow`
These two tools are a combination of fast tables and Parquet. 
`polars` is a DataFrame library built on a fast Rust engine that spreads
work across multiple CPU cores automatically. (much more efficient than one!)
`pyarrow` is the engine that reads/writes .parquet files in python. Together, they
help me loop and cache large tables, so re-parsing doesn't occur more than necessary. 
2. `igraph`
This tool is fast graph math in C: graphs are stored in C arrays under the hood, yet allow
Python interface (shortest paths, connected components, centrality, etc). 
3. `bioframe`
Gives readily made functions to find which genes or peaks overlap using sorted intervals, 
rather than checking every pair by hand. 

**Edge Cases**

As I am working with real data, at times, there are unique case which need to be accounted for. 
Regarding `bf.overlap()`, which looks at the given gene interval and determines overlap, there can be 
unique/edge cases such as: 
1. A gene straddling an anchor boundary (could start inside anchor and extend past end)
**Solution (TBD) - could decide if overlap is fine, or compute an overlap 
2. Spanning across multiple anchor nodes (the same gene would show up as an annotation on multiple nodes)
3. Multiple genes in 1 window
4. Touching intervals (adjacent, zero-length)