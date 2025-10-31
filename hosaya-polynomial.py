
# ==============================================================
#  Description:
#  This SageMath function computes the Hosoya polynomial of the generalized
#  join (joined union) graph R = G[G1, G2, ..., Gn].
# ==============================================================

from sage.all import *

def hosoya_polynomial_generalized_join(G, subgraphs): 
    """
    Computes the Hosoya polynomial of the generalized join graph R = G[G1, G2, ..., Gn],
    according to the formula:

        H(R, x) = M 
                 + [Σ d(G_i, 1) + Σ_{d(v_i,v_j)=1} n_i n_j] x
                 + [Σ (C(n_i,2) - d(G_i,1)) + Σ_{d(v_i,v_j)=2} n_i n_j] x^2
                 + Σ_{k=3}^{D} [Σ_{d(v_i,v_j)=k} n_i n_j] x^k

    where:
        - G  : base graph of order n
        - G1, G2, ..., Gn : subgraphs corresponding to vertices of G
        - d(G_i,1) : number of unordered vertex pairs at distance 1 (edges) within G_i
        - n_i = |V(G_i)| : number of vertices in each subgraph
        - D : diameter of the base graph G
    """

    n = G.order()            # Number of vertices in the base graph
    D = G.diameter()         # Diameter of the base graph
    vertices_G = G.vertices()

    # Total number of vertices in the generalized join
    M = sum(gi.order() for gi in subgraphs)

    # d(G_i,1): number of vertex pairs at distance 1 (edges) in each subgraph
    # In SageMath, .size() returns the number of edges in the graph.
    dGi1 = [gi.size() for gi in subgraphs]

    
    # Coefficient of x (distance 1 pairs)
    
    inter_edges_dist1 = 0
    for i in range(n):
        for j in range(i + 1, n):
            if G.distance(vertices_G[i], vertices_G[j]) == 1:
                inter_edges_dist1 += subgraphs[i].order() * subgraphs[j].order()
    coeff_x = sum(dGi1) + inter_edges_dist1

    
    # Coefficient of x^2 (distance 2 pairs)

    intra_non_edges = sum(binomial(subgraphs[i].order(), 2) - dGi1[i] for i in range(n))

    inter_pairs_dist2 = 0
    for i in range(n):
        for j in range(i + 1, n):
            if G.distance(vertices_G[i], vertices_G[j]) == 2:
                inter_pairs_dist2 += subgraphs[i].order() * subgraphs[j].order()

    # Corrected formula (no extra factor)
    coeff_x2 = intra_non_edges + inter_pairs_dist2

   
    # Coefficients for distances ≥ 3

    poly_dict = {1: coeff_x, 2: coeff_x2}
    for k in range(3, D + 1):
        coeff = 0
        for i in range(n):
            for j in range(i + 1, n):
                if G.distance(vertices_G[i], vertices_G[j]) == k:
                    coeff += subgraphs[i].order() * subgraphs[j].order()
        if coeff > 0:
            poly_dict[k] = coeff

    # Construct the Hosoya polynomial
    R.<x> = PolynomialRing(QQ)
    polynomial = M + sum(c * x**k for k, c in poly_dict.items())

    # Print in readable format
    terms = [f"{M}"] + [f"{c}x^{k}" if k > 1 else f"{c}x" for k, c in poly_dict.items()]
    print("Hosoya polynomial:", " + ".join(terms))

    return polynomial



# Example (for verification and reproducibility)

if __name__ == "__main__":
    G = graphs.PathGraph(4)
    K3 = graphs.CompleteGraph(3)
    P3 = graphs.PathGraph(3)
    C4 = graphs.CycleGraph(4)
    K2 = graphs.CompleteGraph(2)

    H = hosoya_polynomial_generalized_join(G, [K3, P3, C4, K2])
    print("Hosoya polynomial for P4[K3, P3, C4, K2]:")
    print(H)

