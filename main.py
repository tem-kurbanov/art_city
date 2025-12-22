import math
import random
import time
import argparse
import os
import csv
from itertools import combinations
from math import ceil
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# ---------------------------- geometry helpers ---------------------------- #

def euclid(p, q):
    return math.hypot(p[0] - q[0], p[1] - q[1])

def dot(ax, ay, bx, by):
    return ax * bx + ay * by

def angle_ok(p_u, p_v, p_w, min_deg=60.0):
    """
    Check that the angle ∠w u v (between vectors u->w and u->v) is >= min_deg.
    Returns True if OK (i.e., not less than min_deg).
    """
    uxw_x, uxw_y = p_w[0] - p_u[0], p_w[1] - p_u[1]
    uxv_x, uxv_y = p_v[0] - p_u[0], p_v[1] - p_u[1]
    n1 = math.hypot(uxw_x, uxw_y)
    n2 = math.hypot(uxv_x, uxv_y)
    if n1 == 0 or n2 == 0:
        return True  # degenerate; treat as OK
    cos_theta = dot(uxw_x, uxw_y, uxv_x, uxv_y) / (n1 * n2)
    # numerical clamp
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta = math.degrees(math.acos(cos_theta))
    return theta >= min_deg

def _orient(a, b, c):
    """2D orientation (cross product sign)."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

def _on_segment(a, b, c):
    """Assumes a,b,c are colinear; checks if c lies on segment ab."""
    return (min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and
            min(a[1], b[1]) <= c[1] <= max(a[1], b[1]))

def segments_cross_strict(a, b, c, d):
    """
    Proper intersection test: returns True if segments ab and cd intersect
    at an interior point (not just touching at endpoints).
    """
    # Shared endpoints are allowed by the algorithm
    if a == c or a == d or b == c or b == d:
        return False

    o1 = _orient(a, b, c)
    o2 = _orient(a, b, d)
    o3 = _orient(c, d, a)
    o4 = _orient(c, d, b)

    # General case
    if (o1 > 0 and o2 < 0 or o1 < 0 and o2 > 0) and (o3 > 0 and o4 < 0 or o3 < 0 and o4 > 0):
        return True

    # Colinear cases -> treat overlap as crossing (interior overlap only)
    if o1 == 0 and _on_segment(a, b, c):  # c on ab
        return c != a and c != b
    if o2 == 0 and _on_segment(a, b, d):  # d on ab
        return d != a and d != b
    if o3 == 0 and _on_segment(c, d, a):  # a on cd
        return a != c and a != d
    if o4 == 0 and _on_segment(c, d, b):  # b on cd
        return b != c and b != d

    return False

# ---------------------------- graph structure validation ---------------------------- #

def is_strongly_connected(adj_out, n):
    """
    Check if a directed graph is strongly connected (there is a directed path 
    between each pair of vertices).
    
    Args:
        adj_out: dict {i: set of outgoing neighbors} - adjacency list
        n: number of vertices
    
    Returns:
        bool: True if the graph is strongly connected, False otherwise
    """
    if n == 0:
        return True
    if n == 1:
        return True  # Single vertex is trivially strongly connected
    
    # Check if every vertex can reach every other vertex
    # We do this by starting from vertex 0 and checking if it can reach all others,
    # then checking if all others can reach vertex 0 (which implies strong connectivity)
    
    # Forward reachability: can vertex 0 reach all others?
    visited_forward = set()
    queue = [0]
    visited_forward.add(0)
    
    while queue:
        u = queue.pop(0)
        for v in adj_out.get(u, set()):
            if v not in visited_forward:
                visited_forward.add(v)
                queue.append(v)
    
    # If vertex 0 cannot reach all vertices, graph is not strongly connected
    if len(visited_forward) != n:
        return False
    
    # Build reverse adjacency list (incoming edges)
    adj_in_reverse = {i: set() for i in range(n)}
    for u in range(n):
        for v in adj_out.get(u, set()):
            adj_in_reverse[v].add(u)
    
    # Backward reachability: can all others reach vertex 0?
    visited_backward = set()
    queue = [0]
    visited_backward.add(0)
    
    while queue:
        u = queue.pop(0)
        for v in adj_in_reverse.get(u, set()):
            if v not in visited_backward:
                visited_backward.add(v)
                queue.append(v)
    
    # If all vertices cannot reach vertex 0, graph is not strongly connected
    if len(visited_backward) != n:
        return False
    
    return True

# ---------------------------- core generator ---------------------------- #

def generate_graph(n, seed=None, max_tries=1000, min_angle_deg=60.0, 
                   l_correlated=0, m_uncorrelated=0):
    """
    Implements:
      1) Place n nodes uniformly at random in a disk of radius 10*sqrt(n)
      2) Build all potential directed edges with multiple cost criteria:
         - Criterion 0: ceil(Euclidean distance)
         - Criteria 1 to l: correlated with distance (distance + noise)
         - Criteria l+1 to l+m: uncorrelated random values
      3) Sort edges by non-decreasing cost (using first criterion)
      4) Insert a directed edge (i,j) if:
         (a) it doesn't cross any already-inserted edge; and
         (b) at each shared endpoint, it doesn't form an angle < 60° with an existing edge
      5) If any node is isolated (no incoming or outgoing edges), repeat from scratch.
      6) Check that the graph is strongly connected (there is a directed path between each pair of vertices).

    Args:
        n: number of nodes
        seed: random seed
        max_tries: maximum attempts to generate valid graph
        min_angle_deg: minimum angle between edges at a node
        l_correlated: number of criteria correlated with distance
        m_uncorrelated: number of uncorrelated criteria

    Returns:
      positions: dict {i: (x, y)}
      edges: list of (i, j) directed edges
      costs: dict {(i,j): list of costs} where costs[0] is distance, costs[1:l+1] are correlated, costs[l+1:] are uncorrelated
    """
    if seed is not None:
        random.seed(seed)

    R = 10.0 * math.sqrt(n)

    for attempt in range(max_tries):
        # Step 2: positions in disk
        positions = {}
        for i in range(n):
            theta = random.random() * 2 * math.pi
            r = math.sqrt(random.random()) * R  # uniform in disk
            positions[i] = (r * math.cos(theta), r * math.sin(theta))

        # Step 3–5: potential edges with multiple criteria
        P = []
        for i, j in combinations(range(n), 2):
            d = euclid(positions[i], positions[j])
            base_cost = ceil(d)
            
            # Generate multicriteria costs
            costs_list = [base_cost]  # Criterion 0: distance
            
            # Add l correlated criteria (distance + noise)
            for _ in range(l_correlated):
                # Correlated with distance: base_cost + random noise
                noise = random.gauss(0, base_cost * 0.1)  # 10% noise
                correlated_cost = max(1, int(base_cost + noise))
                costs_list.append(correlated_cost)
            
            # Add m uncorrelated criteria (random values)
            for _ in range(m_uncorrelated):
                # Uncorrelated: random value in similar range to distance
                uncorrelated_cost = random.randint(1, max(1, int(base_cost * 2)))
                costs_list.append(uncorrelated_cost)
            
            P.append((costs_list, i, j))
        
        # Sort by first criterion (distance)
        P.sort(key=lambda x: x[0][0])

        E = []                     # list of directed edges (i, j)
        adj_out = {i: set() for i in range(n)}  # outgoing neighbors
        adj_in = {i: set() for i in range(n)}   # incoming neighbors
        costs = {}

        for costs_list, i, j in P:
            pi, pj = positions[i], positions[j]

            # (a) no crossing with existing edges (check both directions)
            crosses = any(segments_cross_strict(pi, pj, positions[u], positions[v]) for (u, v) in E)
            if crosses:
                continue

            # (b) angle constraint at both endpoints: for every neighbor w of i and j
            # Check all neighbors (both incoming and outgoing) for angle constraints
            ok = True
            for w in adj_out[i] | adj_in[i]:  # All neighbors of i
                if not angle_ok(positions[i], positions[j], positions[w], min_deg=min_angle_deg):
                    ok = False
                    break
            if ok:
                for w in adj_out[j] | adj_in[j]:  # All neighbors of j
                    if not angle_ok(positions[j], positions[i], positions[w], min_deg=min_angle_deg):
                        ok = False
                        break
            if not ok:
                continue

            # Insert both directed edges (i, j) and (j, i)
            E.append((i, j))
            E.append((j, i))
            adj_out[i].add(j)
            adj_in[j].add(i)
            adj_out[j].add(i)
            adj_in[i].add(j)
            costs[(i, j)] = costs_list.copy()
            costs[(j, i)] = costs_list.copy()

        # Step 7: accept only if no isolated nodes (must have at least one incoming or outgoing edge)
        if not all(len(adj_out[i]) > 0 or len(adj_in[i]) > 0 for i in range(n)):
            continue
        
        # Step 8: check strong connectivity (directed path between each pair of vertices)
        if not is_strongly_connected(adj_out, n):
            continue
        
        return positions, E, costs

    raise RuntimeError(f"Failed to generate a strongly connected graph without isolated nodes in {max_tries} attempts.")

# ---------------------------- visualization ---------------------------- #

def visualize_graph(positions, edges, costs, title="Generated Directed Graph", show_costs=True, 
                   node_size=100, edge_width=1.5, figsize=(12, 10), criterion_to_show=0):
    """
    Visualize the generated directed graph with nodes, edges, and optionally edge costs.
    
    Args:
        positions: dict {i: (x, y)} - node positions
        edges: list of (i, j) - directed edges
        costs: dict {(i,j): list of costs} - multicriteria edge costs
        title: str - plot title
        show_costs: bool - whether to display edge costs as labels
        node_size: int - size of nodes in the plot
        edge_width: float - width of edges
        figsize: tuple - figure size (width, height)
        criterion_to_show: int - which criterion to display (0=distance, 1=l_correlated, etc.)
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Extract coordinates
    x_coords = [pos[0] for pos in positions.values()]
    y_coords = [pos[1] for pos in positions.values()]
    
    # Draw directed edges with arrows
    for i, j in edges:
        x1, y1 = positions[i]
        x2, y2 = positions[j]
        
        # Calculate arrow position (slightly offset from target node)
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx*dx + dy*dy)
        if length > 0:
            # Normalize and scale arrow offset
            offset = 0.1 * length  # Arrow offset as fraction of edge length
            dx_norm, dy_norm = dx/length, dy/length
            arrow_x = x2 - offset * dx_norm
            arrow_y = y2 - offset * dy_norm
            
            # Draw the edge line
            ax.plot([x1, arrow_x], [y1, arrow_y], 'b-', linewidth=edge_width, alpha=0.7)
            
            # Draw arrowhead
            ax.annotate('', xy=(x2, y2), xytext=(arrow_x, arrow_y),
                       arrowprops=dict(arrowstyle='->', color='blue', lw=edge_width, alpha=0.7))
        
        # Add cost labels if requested
        if show_costs:
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            edge_costs = costs.get((i, j), [0])
            if isinstance(edge_costs, list) and len(edge_costs) > criterion_to_show:
                cost = edge_costs[criterion_to_show]
            else:
                cost = edge_costs[0] if isinstance(edge_costs, list) else edge_costs
            ax.text(mid_x, mid_y, str(cost), ha='center', va='center', 
                   fontsize=8, bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    
    # Draw nodes
    ax.scatter(x_coords, y_coords, c='red', s=node_size, alpha=0.8, edgecolors='black', linewidth=1)
    
    # Add node labels
    for i, (x, y) in positions.items():
        ax.annotate(str(i), (x, y), ha='center', va='center', fontsize=10, fontweight='bold')
    
    # Set equal aspect ratio and add grid
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    # Add title and labels
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('X coordinate', fontsize=12)
    ax.set_ylabel('Y coordinate', fontsize=12)
    
    # Add statistics to the plot
    n_nodes = len(positions)
    n_edges = len(edges)
    
    # Calculate statistics for the displayed criterion
    if n_edges > 0:
        displayed_costs = []
        for edge in edges:
            edge_costs = costs.get(edge, [0])
            if isinstance(edge_costs, list) and len(edge_costs) > criterion_to_show:
                displayed_costs.append(edge_costs[criterion_to_show])
            else:
                displayed_costs.append(edge_costs[0] if isinstance(edge_costs, list) else edge_costs)
        
        total_cost = sum(displayed_costs)
        avg_cost = total_cost / n_edges
        n_criteria = len(list(costs.values())[0]) if costs else 1
    else:
        total_cost = avg_cost = 0
        n_criteria = 1
    
    stats_text = f"Nodes: {n_nodes}, Directed Edges: {n_edges}\nCriteria: {n_criteria}, Showing: {criterion_to_show}\nTotal Cost: {total_cost}, Avg Cost: {avg_cost:.1f}"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout()
    return fig, ax

def save_graph_visualization(positions, edges, costs, filename="graph_visualization.png", **kwargs):
    """
    Save the graph visualization to a file.
    
    Args:
        positions: dict {i: (x, y)} - node positions
        edges: list of (i, j) - edges with i < j
        costs: dict {(i,j): int_cost} - edge costs
        filename: str - output filename
        **kwargs: additional arguments passed to visualize_graph
    """
    fig, ax = visualize_graph(positions, edges, costs, **kwargs)
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Graph visualization saved to {filename}")

def save_graph_csv(positions, edges, costs, filepath, n_criteria=1):
    """
    Save the graph structure to a CSV file.
    
    Format:
    - First line: number of nodes, number of edges
    - Header line: fromId,toId,d,c1,c2,... (d is distance, c1+ are other criteria)
    - Subsequent lines: edge data
    
    Args:
        positions: dict {i: (x, y)} - node positions
        edges: list of (i, j) - directed edges
        costs: dict {(i,j): list of costs} - multicriteria edge costs
        filepath: str - full path to output CSV file
        n_criteria: int - total number of criteria (default: 1)
    """
    n_nodes = len(positions)
    n_edges = len(edges)
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(filepath) if os.path.dirname(filepath) else '.', exist_ok=True)
    
    with open(filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # First line: number of nodes, number of edges
        writer.writerow([n_nodes, n_edges])
        
        # Header line: fromId,toId,d,c1,c2,...
        header = ['fromId', 'toId', 'd']
        for i in range(1, n_criteria):
            header.append(f'c{i}')
        writer.writerow(header)
        
        
        # Write edges sorted by fromId, then toId for consistency
        sorted_edges = sorted(edges)
        for i, j in sorted_edges:
            edge_costs = costs.get((i, j), [0])
            # Ensure we have the right number of costs
            if len(edge_costs) < n_criteria:
                edge_costs = edge_costs + [edge_costs[0]] * (n_criteria - len(edge_costs))
            row = [i, j] + edge_costs[:n_criteria]
            writer.writerow(row)
    
    print(f"Graph saved to {filepath}")

# ---------------------------- cost analysis ---------------------------- #

def analyze_edge_costs(edges, costs, l_correlated=0, m_uncorrelated=0):
    """
    Analyze the distribution of edge costs across all criteria.
    
    Args:
        edges: list of (i, j) directed edges
        costs: dict {(i,j): list of costs} - multicriteria edge costs
        l_correlated: number of correlated criteria
        m_uncorrelated: number of uncorrelated criteria
    
    Returns:
        dict with statistics for each criterion
    """
    n_criteria = 1 + l_correlated + m_uncorrelated
    analysis = {}
    
    for criterion in range(n_criteria):
        criterion_costs = []
        for edge in edges:
            edge_costs = costs.get(edge, [0])
            if len(edge_costs) > criterion:
                criterion_costs.append(edge_costs[criterion])
        
        if criterion_costs:
            criterion_costs = np.array(criterion_costs)
            analysis[criterion] = {
                'values': criterion_costs,
                'mean': np.mean(criterion_costs),
                'std': np.std(criterion_costs),
                'min': np.min(criterion_costs),
                'max': np.max(criterion_costs),
                'range': np.max(criterion_costs) - np.min(criterion_costs),
                'cv': np.std(criterion_costs) / np.mean(criterion_costs) if np.mean(criterion_costs) > 0 else 0,
                'type': 'distance' if criterion == 0 else 
                       'correlated' if 1 <= criterion <= l_correlated else 'uncorrelated'
            }
    
    return analysis

def visualize_cost_distributions(edges, costs, l_correlated=0, m_uncorrelated=0, 
                                figsize=(15, 10), bins=20):
    """
    Create visualizations showing the distribution of edge costs across all criteria.
    
    Args:
        edges: list of (i, j) directed edges
        costs: dict {(i,j): list of costs} - multicriteria edge costs
        l_correlated: number of correlated criteria
        m_uncorrelated: number of uncorrelated criteria
        figsize: figure size
        bins: number of bins for histograms
    """
    analysis = analyze_edge_costs(edges, costs, l_correlated, m_uncorrelated)
    n_criteria = len(analysis)
    
    if n_criteria == 0:
        print("No criteria to analyze")
        return None, None
    
    # Create subplots
    fig, axes = plt.subplots(2, n_criteria, figsize=figsize)
    if n_criteria == 1:
        axes = axes.reshape(2, 1)
    
    fig.suptitle('Edge Cost Distributions Across All Criteria', fontsize=16, fontweight='bold')
    
    for criterion, stats in analysis.items():
        col = criterion
        
        # Histogram
        ax_hist = axes[0, col]
        ax_hist.hist(stats['values'], bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
        ax_hist.set_title(f'Criterion {criterion} ({stats["type"]})\n'
                         f'Mean: {stats["mean"]:.2f}, Std: {stats["std"]:.2f}, CV: {stats["cv"]:.3f}')
        ax_hist.set_xlabel('Cost Value')
        ax_hist.set_ylabel('Frequency')
        ax_hist.grid(True, alpha=0.3)
        
        # Box plot
        ax_box = axes[1, col]
        box_data = ax_box.boxplot(stats['values'], patch_artist=True)
        box_data['boxes'][0].set_facecolor('lightcoral')
        ax_box.set_title(f'Criterion {criterion} Box Plot\n'
                        f'Range: {stats["min"]:.0f} - {stats["max"]:.0f}')
        ax_box.set_ylabel('Cost Value')
        ax_box.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig, axes

def print_cost_analysis(edges, costs, l_correlated=0, m_uncorrelated=0):
    """
    Print detailed analysis of edge cost distributions.
    """
    analysis = analyze_edge_costs(edges, costs, l_correlated, m_uncorrelated)
    
    print("\n" + "="*80)
    print("EDGE COST DISTRIBUTION ANALYSIS")
    print("="*80)
    
    for criterion, stats in analysis.items():
        print(f"\nCriterion {criterion} ({stats['type'].upper()}):")
        print(f"  Count: {len(stats['values'])} edges")
        print(f"  Mean:  {stats['mean']:.2f}")
        print(f"  Std:   {stats['std']:.2f}")
        print(f"  Min:   {stats['min']:.0f}")
        print(f"  Max:   {stats['max']:.0f}")
        print(f"  Range: {stats['range']:.0f}")
        print(f"  CV:    {stats['cv']:.3f} (coefficient of variation)")
        
        # Uniformity assessment
        if stats['cv'] < 0.1:
            uniformity = "Very uniform"
        elif stats['cv'] < 0.2:
            uniformity = "Fairly uniform"
        elif stats['cv'] < 0.5:
            uniformity = "Moderately varied"
        else:
            uniformity = "Highly varied"
        print(f"  Uniformity: {uniformity}")
    
    # Correlation analysis
    if l_correlated > 0:
        print(f"\nCORRELATION ANALYSIS:")
        distance_values = analysis[0]['values']
        for criterion in range(1, l_correlated + 1):
            if criterion in analysis:
                correlated_values = analysis[criterion]['values']
                correlation = np.corrcoef(distance_values, correlated_values)[0, 1]
                print(f"  Criterion {criterion} vs Distance: {correlation:.3f}")

# ---------------------------- example usage ---------------------------- #

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate and visualize VRP graphs')
    parser.add_argument('-n', '--nodes', type=int, default=50, 
                       help='Number of nodes in the graph (default: 50)')
    parser.add_argument('--no-show', action='store_true', 
                       help='Do not display the plot window (useful for headless environments)')
    parser.add_argument('--no-save', action='store_true', 
                       help='Do not save the visualization to file')
    parser.add_argument('--no-costs', action='store_true', 
                       help='Do not show edge costs in the visualization')
    parser.add_argument('-l', '--correlated', type=int, default=0,
                       help='Number of criteria correlated with distance (default: 0)')
    parser.add_argument('-m', '--uncorrelated', type=int, default=0,
                       help='Number of uncorrelated criteria (default: 0)')
    parser.add_argument('--show-criterion', type=int, default=0,
                       help='Which criterion to display in visualization (0=distance, 1=first correlated, etc.)')
    parser.add_argument('--analyze-costs', action='store_true',
                       help='Show detailed cost distribution analysis')
    parser.add_argument('--cost-plots', action='store_true',
                       help='Generate cost distribution plots')
    parser.add_argument('-o', '--output-dir', type=str, default=None,
                       help='Directory to save the graph CSV file (default: current directory)')
    parser.add_argument('-g', '--num-graphs', type=int, default=1,
                       help='Number of graphs to generate (default: 1)')
    
    args = parser.parse_args()
    
    n = args.nodes
    if n < 2:
        print("Error: Number of nodes must be at least 2")
        exit(1)
    
    if args.num_graphs < 1:
        print("Error: Number of graphs must be at least 1")
        exit(1)
    
    n_criteria = 1 + args.correlated + args.uncorrelated
    
    # Base seed from current time
    base_seed = int(time.time() * 1000) % 100000  # Use milliseconds, modulo for reasonable size
    
    # Generate multiple graphs
    for graph_idx in range(args.num_graphs):
        print(f"\n{'='*80}")
        print(f"Generating graph {graph_idx + 1} of {args.num_graphs}")
        print(f"{'='*80}")
        
        # Generate unique seed for each graph
        graph_seed = (base_seed + graph_idx) % 100000
        
        positions, edges, costs = generate_graph(n, seed=graph_seed, 
                                               l_correlated=args.correlated, 
                                               m_uncorrelated=args.uncorrelated)
        
        print(f"n={n}, |E|={len(edges)}, criteria={n_criteria}, seed={graph_seed}")
        print(f"  - Distance (criterion 0)")
        if args.correlated > 0:
            print(f"  - {args.correlated} correlated criteria (1-{args.correlated})")
        if args.uncorrelated > 0:
            print(f"  - {args.uncorrelated} uncorrelated criteria ({args.correlated+1}-{n_criteria-1})")
        
        # Print a few edges with costs
        for k, (i, j) in enumerate(edges[:10], 1):
            edge_costs = costs[(i,j)]
            cost_str = f"[{', '.join(map(str, edge_costs))}]" if len(edge_costs) > 1 else str(edge_costs[0])
            print(f"{k:2d}. ({i}, {j})  costs={cost_str}")
        
        # Save graph to CSV if output directory is specified
        if args.output_dir is not None:
            # Create informative filename with graph index if generating multiple
            if args.num_graphs > 1:
                graph_filename = f"vrp_multicriteria_n{n}_l{args.correlated}_m{args.uncorrelated}_seed{graph_seed}_graph{graph_idx+1:04d}.csv"
            else:
                graph_filename = f"vrp_multicriteria_n{n}_l{args.correlated}_m{args.uncorrelated}_seed{graph_seed}.csv"
            graph_filepath = os.path.join(args.output_dir, graph_filename)
            save_graph_csv(positions, edges, costs, graph_filepath, n_criteria)
        
        # Visualize the graph (only for first graph or if explicitly requested)
        if graph_idx == 0 or args.num_graphs == 1:
            print("\nGenerating graph visualization...")
            fig, ax = visualize_graph(positions, edges, costs, 
                                    title=f"VRP Multicriteria Graph (n={n}, |E|={len(edges)}, criteria={n_criteria}, seed={graph_seed})", 
                                    show_costs=not args.no_costs,
                                    criterion_to_show=args.show_criterion)
            
            # Save the visualization unless disabled
            if not args.no_save:
                if args.num_graphs > 1:
                    vis_filename = f"vrp_multicriteria_n{n}_l{args.correlated}_m{args.uncorrelated}_seed{graph_seed}_graph{graph_idx+1:04d}.png"
                else:
                    vis_filename = f"vrp_multicriteria_n{n}_l{args.correlated}_m{args.uncorrelated}_seed{graph_seed}.png"
                save_graph_visualization(positions, edges, costs, 
                                       filename=vis_filename,
                                       title=f"VRP Multicriteria Graph (n={n}, |E|={len(edges)}, criteria={n_criteria}, seed={graph_seed})",
                                       criterion_to_show=args.show_criterion)
            
            # Show the plot unless disabled (only for first graph)
            if not args.no_show and graph_idx == 0:
                plt.show()
            else:
                plt.close(fig)
        
        # Cost analysis (only for first graph or if explicitly requested)
        if args.analyze_costs and (graph_idx == 0 or args.num_graphs == 1):
            print_cost_analysis(edges, costs, args.correlated, args.uncorrelated)
        
        # Cost distribution plots (only for first graph or if explicitly requested)
        if args.cost_plots and (graph_idx == 0 or args.num_graphs == 1):
            print("\nGenerating cost distribution plots...")
            fig_costs, axes_costs = visualize_cost_distributions(edges, costs, args.correlated, args.uncorrelated)
            
            if fig_costs is not None:
                # Save cost distribution plots
                if args.num_graphs > 1:
                    cost_plot_filename = f"cost_distributions_n{n}_l{args.correlated}_m{args.uncorrelated}_seed{graph_seed}_graph{graph_idx+1:04d}.png"
                else:
                    cost_plot_filename = f"cost_distributions_n{n}_l{args.correlated}_m{args.uncorrelated}_seed{graph_seed}.png"
                fig_costs.savefig(cost_plot_filename, dpi=300, bbox_inches='tight')
                print(f"Cost distribution plots saved to {cost_plot_filename}")
                
                if not args.no_show and graph_idx == 0:
                    plt.show()
                else:
                    plt.close(fig_costs)
    
    print(f"\n{'='*80}")
    print(f"Successfully generated {args.num_graphs} graph(s)")
    print(f"{'='*80}")
