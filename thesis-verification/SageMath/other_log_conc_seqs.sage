import sys
load("bivariate_chromatic.sage")

# Verification used for section 6.2 of my master's thesis (titled "Terms in the Deletion-Contraction Recurrence")

# Lower and upper bounds on number of vertices
NUM_VERTICES_LOWER = 1
NUM_VERTICES_UPPER = 8


def add_seqs(s1,s2):
    """ Returns the component-wise addition of s1 and s2. """
    if len(s1) != len(s2):
        print "add_seqs: s1,s2 have different length"
        return None
    return [s1[j] + s2[j] for j in range(len(s1))]

def sub_seqs(s1,s2):
    """ Returns the component-wise subtraction of s1 and s2. """
    if len(s1) != len(s2):
        print "sub_seqs: s1,s2 have different length"
        return None
    return [s1[j] - s2[j] for j in range(len(s1))]

def is_log_concave(seq):
    """ Returns true if the sequence seq is log-concave, false otherwise. """
    for j in range(1,len(seq)-1):
        if seq[j]^2 < seq[j-1]*seq[j+1]:
            return false
    return true


def test_graph(g6str,kmax):
    """ Compute the sequences and test log-concavity """
    G = Graph(g6str)
    P = G.bivariate_chrom_poly()
    for e in G.edges():
        Gdel = G.copy()
        Gcon = G.copy()
        Gext = G.copy()

        Gdel.delete_edge(e)
        Gcon.contract_edge(e)
        Gext.delete_vertices([e[0],e[1]])

        Pdel = Gdel.bivariate_chrom_poly()
        Pcon = Gcon.bivariate_chrom_poly()
        Pext = Gext.bivariate_chrom_poly()

        for k in range(1,kmax):
            S = [P(x=k,y=j) for j in range(k+1)]
            Sdel = [Pdel(x=k,y=j) for j in range(k+1)]
            Scon = [Pcon(x=k,y=j) for j in range(k+1)]
            Sext = [(k-j)*Pext(x=k,y=j) for j in range(k+1)]

            if not is_log_concave(sub_seqs(Sdel,Scon)):
                print "Sdel - Scon is not log-concave for" + g6str
                return
            if not is_log_concave(sub_seqs(Sext, Scon)):
                print "-Scon + Sext is not log-concave for" + g6str
                return
            if not is_log_concave(add_seqs(Sdel,Sext)):
                print "Sdel + Sext is not log-concave for" + g6str
                return
    # print("Done " + g6str) # Debugging


for n in range(NUM_VERTICES_LOWER, NUM_VERTICES_UPPER):
    with open("../graph_data/connected/graphs_{0}.g6".format(n), 'r') as f:
        count = 0
        print "Starting graphs on {0} vertices".format(n)
        for g6str in f:
            if (count % 100 == 0):
                print "{0:5} |".format(count),
            g6str = g6str.strip()
            test_graph(g6str, n+1)
            count += 1
            if (count % 10 == 0):
                print "-",
                sys.stdout.flush()
                if (count % 100 == 0):
                    print "|"
    print ""
    print "Finished graphs on {0} vertices".format(n)
