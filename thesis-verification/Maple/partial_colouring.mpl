interface(quiet=true):

with(combinat,powerset):
with(GraphTheory):


### How to use this file:

## A file should start with something like:

# interface(quiet=true):  # optional
# read("partial_colouring.mpl")


## To iterate over graphs:

# iter := ConnectedGraphIter(n)
# do:
#     G := iter();
#     if G = FAIL then
#         break;
#     end if;

#     # Do stuff

# end do:




# Returns exp(sum(t_v, v in V(G))) + x*sum(t^S, S stable set in G) if groundSet is V(G)
Qseries := proc(groundSet::set, stableSets::set, $)
    return exp(add(t[v], v in groundSet)) + x*add(mul(t[v], v in S), S in stableSets);
end proc:


# Returns true if no element of E is a subset of S (i.e. S is stable)
isStable := proc(S::set, E::set(set), $)::boolean;
    return not ormap((s) -> s subset S, E);
end proc:


# Some lists to avoid repeating common computations
_vSet[0] := {}:    # List of vertex sets
_vPowSet[0] := {}: # List of powersets of vertex sets

# Initialize a graph G for future use of the Bkpoly function
# Note: Using the "Z series" to compute B_k is the fastest way I found for computing B_k
#       In particular, computing the bivariate chromatic polynomial can be quite slow
InitGraph := proc(G::Graph, $)
    global _vSet, _vPowSet;
    local n, stableSets, Z;

    n := NumberOfVertices(G);
    if not type(_vSet[n], set) then
        _vSet[n] := {seq(i, i=1..n)};
        _vPowSet[n] := powerset(_vSet[n]);
    end if;

    stableSets := select(isStable, _vPowSet[n], Edges(G));
    Z := 1 / (1 - y*Qseries(_vSet[n], stableSets));
    SetGraphAttribute(G, "Zseries"=Z);
end proc:


# Compute the Bivariate Chromatic Polynomial of G
BivariateChromaticPoly := proc(G::Graph, x::name, y::name, $)::polynom(integer);
    global _vSet, _vPowSet;
    local P, n, W;
    n := NumberOfVertices(G);
    P := 0;

    # Take care of this case
    if (n = 0) then
        return 1;
    elif (n = 1) then
        return x;
    end if;

    if not type(_vSet[n], set) then
        _vSet[n] := {seq(i, i=1..n)};
        _vPowSet[n] := powerset(_vSet[n]);
    end if;

    P := add((x-y)^(n - nops(W))*ChromaticPolynomial(InducedSubgraph(G,W),y), W in _vPowSet[n]);
    return simplify(P)
end proc:


# Returns the coefficient of the the monomial which is the product of the elements of the set S raised to the given powers
coeffList := proc(expr, S::set(list), $)
    return foldl((Z, x) -> coeftayl(Z, x[1] = 0, x[2]), expr, op(S));
end proc:


# Get the B_k polynomial for the graph G (assumes that InitGraph has already been called on G)
Bkpoly := proc(G::Graph, k::nonnegint, $)
    local Z;
    Z := GetGraphAttribute(G, "Zseries");
    return simplify(coeffList(coeftayl(Z, y = 0, k), {seq([t[v], 1], v in _vSet[NumberOfVertices(G)])}));
end proc:


# Compute the B_k polynomial using the bivariate chromatic polynomial (used for testing other procedures)
BpolyFromBiChromPoly := proc(P::polynom, kval::nonnegint, $)::polynom;
    return add(binomial(kval,j) * eval(P,[x=kval,y=j]) * x^j, j=0..kval);
end proc:



### These returns iterators for the graphs of a desired type
### The iterator will return the next graph or FAIL

# For connected graphs on n vertices
# (The graph will have an attribute named "g6string" listing
#  the g6 representation it came from)
ConnectedGraphIter := proc(n::nonnegint, $)::string;
    local path,fd;
    path := cat("/Users/glmckay/Documents/Research/graph_data/connected/graphs_", convert(n, string), ".g6");
    # return ImportGraph(path, "graph6", output=iterator);
    fd := fopen(path, READ, BINARY);
    return proc()
        local g6str,G;
        if (fd < 0) then
            return FAIL;
        end if;
        g6str := readline(fd);
        if (g6str = 0) then
            fclose(fd);
            fd := -1;
            return FAIL;
        end if;
        G := ConvertGraph(g6str);
        SetGraphAttribute(G, "g6string"=g6str);
        InitGraph(G);
        return G;
    end proc;
end proc:

# For connected graphs on n vertices whose complement is triangle-free
# (for now this does not store the g6string as an attribute or call InitGraph)
TriFreeComplGraphIter := proc(n::nonnegint, $)::string;
    local path;
    path := cat("/Users/glmckay/Documents/Research/graph_data/tri_free_compl/graphs_", convert(n, string), ".g6");
    return ImportGraph(path, "graph6", output=iterator);
end proc:

# For Biconnected (2-connected or K_2) graphs on n vertices
# (for now this does not store the g6string as an attribute or call InitGraph)
BiconnectedGraphIter := proc(n::nonnegint, $)::string;
    local path;
    path := cat("/Users/glmckay/Documents/Research/graph_data/biconnected/graphs_", convert(n, string), ".g6");
    return ImportGraph(path, "graph6", output=iterator);
end proc:


# Version of the 'DeleteEdge' function that is safe for use with Bkpoly
# (uses the 'inplace=false' option)
SafeDeleteEdge:= proc(G::Graph, e, $)::Graph;
    local H;
    H := DeleteEdge(G, e, inplace=false);
    InitGraph(H);
    return H;
end proc;


# Version of the 'Contract' function that is safe for use with Bkpoly
SafeContract := proc(G::Graph, e, $)::Graph;
    local n,H;
    n := NumberOfVertices(G) - 1;
    H := RelabelVertices(Contract(G,e), [seq(i, i=1..n)]);
    InitGraph(H);
    return H;
end proc;


# Version of the 'Contract' function that is safe for use with Bkpoly
SafeDeleteVertex := proc(G::Graph, v, $)::Graph;
    local n,H;
    n := NumberOfVertices(G) - nops(v);
    if (n = 0) then
        H := Graph();
    else
        H := RelabelVertices(DeleteVertex(G,v), [seq(i, i=1..n)]);
    end if;
    InitGraph(H);
    return H;
end proc;
