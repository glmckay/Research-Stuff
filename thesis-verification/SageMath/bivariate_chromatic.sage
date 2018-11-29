# Add some functions to SageMath's graph class
# I think in the most up-to-date version of SageMath, a contract_edge function implemented by default

y = var('y')


def contract_edge(G,e):
    u = e[0]
    v = e[1]
    nbrs = set(G.neighbors(u)).union(G.neighbors(v))-set([u,v])
    G.delete_vertices([u,v])
    x = G.add_vertex()
    G.add_edges([(x,y) for y in nbrs])


def bivariate_chrom_poly(G):
    p = 0
    n = len(G.vertices())
    for W in Subsets(G.vertices()):
        p = p + (x-y)^(n - len(W)) * G.subgraph(W).chromatic_polynomial().subs(x=y)
    return p.expand().simplify()


def B_k(G,k):
    p = bivariate_chrom_poly(G)
    s = 0
    for j in range(k+1):
        s += binomial(k,j)*p.subs(x=k,y=j)*y^j
    return s.expand().simplify()


Graph.contract_edge = contract_edge
Graph.bivariate_chrom_poly = bivariate_chrom_poly
Graph.B_k = B_k
