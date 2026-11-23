import os
import bioframe as bf
import igraph as ig
import numpy as np
import pandas as pd
import polars as pl

LOOPS_TXT = "LI_fithic_obs-count_full_Q05.tx"
LOOPS_PARQUET = "fithic_loops.parquet"
GENES_BED = "GRCh38-cCREs.ELS.bed"  # Replace with your gene BED/GTF file
TARGET_CHROMS = ["chr21"]  # Change to 'None' to run ALL chromosomes at once

# 1. Convert raw .tx to binary Parquet once
if not os.path.exists(LOOPS_PARQUET):
  (
      pl.scan_csv(LOOPS_TXT, separator="\t")
      .select([
          "rchrom",
          "rstart",
          "rend",
          "rwindow",
          "rtad",
          "cchrom",
          "cstart",
          "cend",
          "cwindow",
          "ctad",
          "obs_count",
          "q",
      ])
      .sink_parquet(LOOPS_PARQUET)
      # executes into chunks, writes result directly to .parquet w/o loading entire dataset
  )

# 2. Lazy-query Parquet (only reads filtered rows/columns into RAM)
q_loops = pl.scan_parquet(LOOPS_PARQUET).filter(pl.col("q") <= 0.05)

if TARGET_CHROMS is not None:
  q_loops = q_loops.filter(
      pl.col("rchrom").is_in(TARGET_CHROMS)
      & pl.col("cchrom").is_in(TARGET_CHROMS)
  )

loops_pl = q_loops.collect()

# 3. Build Nodes DataFrame in Polars
r_nodes = loops_pl.select([
    pl.col("rwindow").alias("node"),
    pl.col("rchrom").alias("chrom"),
    pl.col("rstart").alias("start"),
    pl.col("rend").alias("end"),
    pl.col("rtad").alias("tad"),
])
c_nodes = loops_pl.select([
    pl.col("cwindow").alias("node"),
    pl.col("cchrom").alias("chrom"),
    pl.col("cstart").alias("start"),
    pl.col("cend").alias("end"),
    pl.col("ctad").alias("tad"),
])

nodes_pd = (
    pl.concat([r_nodes, c_nodes])
    .unique(subset=["node"])
    .sort(["chrom", "start"])
    .to_pandas()
)

# 4. Overlap ALL genes onto the anchor nodes
# Assuming a 6-column BED: chrom, start, end, gene_name, score, strand
genes_pd = pd.read_csv(
    GENES_BED,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "gene_name", "score", "strand"],
    usecols=["chrom", "start", "end", "gene_name", "strand"],
)
if TARGET_CHROMS is not None:
  genes_pd = genes_pd[genes_pd["chrom"].isin(TARGET_CHROMS)]

overlaps = bf.overlap(nodes_pd, genes_pd, how="left", suffixes=("", "_gene"))

gene_agg = (
    overlaps.dropna(subset=["gene_name_gene"])
    .groupby("node", as_index=False)
    .agg(
        gene_count=("gene_name_gene", "nunique"),
        genes=("gene_name_gene", lambda s: ",".join(sorted(set(s)))),
    )
)

nodes_pd = nodes_pd.merge(gene_agg, on="node", how="left")
nodes_pd["gene_count"] = nodes_pd["gene_count"].fillna(0).astype(int)
nodes_pd["genes"] = nodes_pd["genes"].fillna("")

# 5. Build Edges DataFrame & C-backed igraph.Graph
edges_pd = (
    loops_pl.group_by(["rwindow", "cwindow"])
    .agg([
        pl.col("obs_count").sum().alias("obs_count"),
        pl.col("q").min().alias("min_q"),
    ])
    .to_pandas()
)
edges_pd["neg_log10_q"] = -np.log10(edges_pd["min_q"].clip(lower=1e-300))

# igraph.Graph.DataFrame builds the entire C graph directly from DataFrames!
g = ig.Graph.DataFrame(
    edges=edges_pd,
    directed=False,
    vertices=nodes_pd,  # First column ('node') automatically becomes vertex 'name'
    use_vids=False,
)

print(
    f"Built igraph with {g.vcount():,} vertices and {g.ecount():,} edges in C!"
)
# Example query: inspect genes on the first few vertices
for v in g.vs[:5]:
  if v["gene_count"] > 0:
    print(f"Anchor {v['name']} ({v['chrom']}:{v['start']}-{v['end']}) -> Genes: {v['genes']}")