# VRP Multicriteria Graph Generator

A Python tool for generating strongly connected directed graphs for Vehicle Routing Problem (VRP) research with support for multiple cost criteria.

## Overview

This tool generates geometric directed graphs suitable for VRP instances with the following characteristics:

- **Strongly Connected**: There exists a directed path between every pair of vertices
- **Non-crossing Edges**: Edges do not intersect at interior points
- **Angle Constraints**: Edges meeting at a node form angles of at least 60°
- **Multiple Cost Criteria**: Supports distance-based and additional correlated/uncorrelated cost criteria
- **Geometric Layout**: Nodes are placed uniformly at random in a disk

## Features

- Generate strongly connected directed graphs with customizable node counts
- Support for multiple cost criteria (distance + correlated + uncorrelated)
- Automatic validation of graph structure (strong connectivity, no isolated nodes)
- Visualization of generated graphs with matplotlib
- Export graphs to CSV format for use in optimization algorithms
- Batch generation of multiple graphs
- Cost distribution analysis and visualization

## Requirements

- Python 3.6+
- numpy
- matplotlib

## Installation

1. Clone or download this repository
2. Install required dependencies:

```bash
pip install numpy matplotlib
```

## Usage

### Basic Usage

Generate a graph with default settings (50 nodes):

```bash
python main.py
```

### Generate a Graph with Specific Parameters

```bash
python main.py -n 100 -l 2 -m 1 -o ./graphs
```

This generates a graph with:
- `-n 100`: 100 nodes
- `-l 2`: 2 correlated cost criteria
- `-m 1`: 1 uncorrelated cost criterion
- `-o ./graphs`: Save CSV to `./graphs` directory

### Generate Multiple Graphs

```bash
python main.py -n 100 -g 10 -o ./graphs
```

Generates 10 graphs with 100 nodes each, saving all CSV files to `./graphs`.

### Headless Mode (No Visualization)

```bash
python main.py -n 200 -g 50 -o ./graphs --no-show --no-save
```

Useful for batch generation without visualization overhead.

### Cost Analysis

```bash
python main.py -n 100 --analyze-costs --cost-plots
```

Shows detailed cost distribution statistics and generates distribution plots.

## Command-Line Arguments

| Argument | Short | Type | Default | Description |
|----------|-------|------|---------|-------------|
| `--nodes` | `-n` | int | 50 | Number of nodes in the graph |
| `--correlated` | `-l` | int | 0 | Number of criteria correlated with distance |
| `--uncorrelated` | `-m` | int | 0 | Number of uncorrelated criteria |
| `--output-dir` | `-o` | str | None | Directory to save graph CSV files |
| `--num-graphs` | `-g` | int | 1 | Number of graphs to generate |
| `--show-criterion` | | int | 0 | Which criterion to display (0=distance, 1+=others) |
| `--no-show` | | flag | False | Do not display plot window |
| `--no-save` | | flag | False | Do not save visualization PNG files |
| `--no-costs` | | flag | False | Do not show edge costs in visualization |
| `--analyze-costs` | | flag | False | Show detailed cost distribution analysis |
| `--cost-plots` | | flag | False | Generate cost distribution plots |

## Output Format

### CSV File Format

Graphs are saved in CSV format with the following structure:

```
n_nodes,n_edges

fromId,toId,d,c1,c2,...
i,j,cost_0,cost_1,cost_2,...
...
```

- **First line**: Number of nodes, number of edges
- **Blank line**
- **Header line**: Column names (`fromId`, `toId`, `d` for distance, `c1`, `c2`, ... for additional criteria)
- **Blank line**
- **Data rows**: Each row represents a directed edge with its costs

Example:
```
100,276

fromId,toId,d,c1
19,85,1,1
85,19,1,1
12,39,2,1
39,12,2,1
...
```

### File Naming

- **Single graph**: `vrp_multicriteria_n{node_count}_l{correlated}_m{uncorrelated}_seed{seed}.csv`
- **Multiple graphs**: `vrp_multicriteria_n{node_count}_l{correlated}_m{uncorrelated}_seed{seed}_graph{index:04d}.csv`

## Algorithm Details

The graph generation algorithm follows these steps:

1. **Node Placement**: Place `n` nodes uniformly at random in a disk of radius `10√n`

2. **Edge Cost Calculation**: For each potential edge `(i,j)`:
   - **Criterion 0 (Distance)**: `ceil(Euclidean distance)`
   - **Correlated Criteria (1 to l)**: Distance + Gaussian noise (10% variance)
   - **Uncorrelated Criteria (l+1 to l+m)**: Random values in similar range

3. **Edge Selection**: Sort potential edges by distance (criterion 0), then add edges that:
   - Do not cross any existing edge (strict intersection test)
   - Satisfy angle constraints (≥60° at each endpoint)
   - Are added in both directions `(i,j)` and `(j,i)`

4. **Validation**: Accept graph only if:
   - No isolated nodes (each node has at least one incoming or outgoing edge)
   - Graph is strongly connected (directed path exists between every pair of vertices)

5. **Retry**: If validation fails, restart from step 1 (up to 1000 attempts)

## Graph Properties

- **Strong Connectivity**: Ensures every vertex is reachable from every other vertex
- **Bidirectional Edges**: Each edge `(i,j)` has a corresponding edge `(j,i)` with the same costs
- **Planar-like Structure**: No edge crossings (except at shared endpoints)
- **Geometric Constraints**: Minimum 60° angle between edges at each node

## Examples

### Example 1: Simple Distance-Only Graph

```bash
python main.py -n 50 -o ./output
```

Generates a graph with 50 nodes, only distance costs, saves to `./output`.

### Example 2: Multi-Criteria Graph

```bash
python main.py -n 100 -l 2 -m 1 -o ./output
```

Generates a graph with:
- 100 nodes
- Distance (criterion 0)
- 2 correlated criteria (criteria 1-2)
- 1 uncorrelated criterion (criterion 3)

### Example 3: Batch Generation

```bash
python main.py -n 200 -g 20 -o ./batch_output --no-show
```

Generates 20 graphs with 200 nodes each, no visualization windows.

### Example 4: With Analysis

```bash
python main.py -n 100 -l 1 -m 1 --analyze-costs --cost-plots
```

Generates a graph and shows detailed cost statistics and distribution plots.

## Visualization

The tool generates visualizations showing:
- Node positions and labels
- Directed edges with arrows
- Edge costs (optional)
- Graph statistics (node count, edge count, criteria info)

Visualizations are saved as PNG files unless `--no-save` is used.

## Notes

- Graph generation may take longer for larger node counts or when many criteria are specified
- The algorithm may need multiple attempts to generate a valid strongly connected graph
- Seeds are automatically generated based on current time, but can be made deterministic by modifying the code

## License

This project is provided as-is for research and educational purposes.

