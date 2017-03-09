#! /usr/bin/env python3

import argparse
import hashlib
import random
import sys


class VertexNames():
    def __init__(self):
        self._names = dict()
        self._ids = dict()
        self._current_id = 0    # replace long string node names in input with unique identifiers

    def get_id(self, name):
        if name in self._names:
            return self._names[name]
        else:
            self._names[name] = self._current_id
            id = self._current_id
            self._ids[id] = name
            self._current_id += 1
            return id

    def get_name(self, id):
        return self._ids[id]


class Graph():
    def __init__(self):
        self._adjacencies = dict()

    def __repr__(self):
        return '<Graph ' + str(self._adjacencies) + '>'

    def add_edge(self, a, b):
        self.add_directed_edge(a, b)
        self.add_directed_edge(b, a)

    def add_directed_edge(self, a, b):
        if a in self._adjacencies:
            self._adjacencies[a].add(b)
        else:
            self._adjacencies[a] = set([b])

    def vertices(self):
        return self._adjacencies.keys()

    def n_vertices(self):
        return len(self.vertices())

    def adjacencies(self, v):
        return self._adjacencies[v]

    def degree(self, v):
        return len(self.adjacencies(v))

    def remove_high_degree_nodes(self, max_degree):
        bad_ones = set()
        poor_neighbors = set()  # Keep track of neighbors to "bad ones"
        
        # First find the nodes to delete
        for v in self.vertices():
            if len(self.adjacencies(v)) > max_degree:
                bad_ones.add(v)
                poor_neighbors |= self.adjacencies(v) # Add v's neighboors to the poor neighbors

        # A bad vertex cannot be a poor neighbor
        poor_neighbors -= bad_ones

        # Then delete bad vertices
        for v in bad_ones:
            del self._adjacencies[v]

        # All high-degree nodes are now removed, but there are still edges pointing at these,
        # so we find these edges and remove them
        for u in poor_neighbors:
            self._adjacencies[u] -= bad_ones
        
        return len(bad_ones)


    def remove_singletons(self):
        singletons = []
        for v in self._adjacencies:
            if len(self._adjacencies[v]) == 0:
                singletons.append(v)
        for v in singletons:
            del self._adjacencies[v]
        return len(singletons)
    
    def component(self, v):
        '''
        Return set of nodes reachable from v. 
        Based on code found on the net.
        '''
        c = set()               # The component we are "collecting"
        nodes = set([v])        # Current set of collected nodes
        while nodes:
            v = nodes.pop()     # Grab one of the nodes
            c.add(v)            # Add it to output collection
            nodes |= self._adjacencies[v] - c # Add neighbors that have not been seen before to the current set
        return c

    def connected_components(self):
        '''
        ?
        '''
        seen = set()
        for node in self._adjacencies:
            if node not in seen:
                c = self.component(node)
                seen |= c
                yield c


    def component_limited(self, v, maxdegree):
        '''
        Return set of nodes reachable from v. Restrict to components of limited degree, see
        degree_limited_connected_components()!

        v is expected to have low degree.

        Based on code found on the net.
        '''
        c = set()               # The component we are "collecting"
        nodes = set([v])        # Current set of collected nodes
        while nodes:
            v = nodes.pop()     # Grab one of the nodes
            c.add(v)            # Add it to output collection
            if self.degree(v) <= maxdegree:
                nodes |= self._adjacencies[v] - c # Add neighbors that have not been seen before to the current set
            # Note that a high degree node is added to the output, but not for further investigation!
        return c


    def degree_limited_connected_components(self, maxdegree):
        '''
        Compute connected components, but stop extending components when a high-degree node is found.
        The idea is that we want components that represent contigs, but do not include repeat structures.
        Repeats create high-degree clusters. However, we want to be able to create contigs that extend 
        into repeats, so a "last" node with many neighbors is OK. This would correspond to a node/contig
        that has only a part of it in a repeat structure (many neighbors), but another part with 
        non-repeat structure (some neighbors have low degree).
        '''
        seen = set()
        for node in self._adjacencies:
            if node not in seen and self.degree(node) <= maxdegree:
                c = self.component_limited(node, maxdegree)
                seen |= c
                yield c


    def output_component_dot(self, g, names, filename, c):
        '''
        Write component c to filename in dot format.
        '''
        with open(filename, "w") as oh:
            oh.write("graph random_component_size_" + str(len(c)) + " {\n")
            for node in c:
                for a in g.adjacencies(node):
                    if node < a:
                        node_name = names.get_name(node)
                        a_name = names.get_name(a)
                        oh.write("  " + str(node) + " -- " + str(a) + "; // " + node_name + ", " + a_name + "\n")
            oh.write("}\n")




def read_graph(f, quiet):
    '''
    Reads a graph in "M4" format, used by dfp_overlap.
    '''
    names = VertexNames()
    g = Graph()
    ignore_vertices = dict()             # Store ids of contigs contained within other contigs
    candidate_edges = []                 # Edges that we might want to have. Will have to check the ignore_vertices before committing.

    with open(f, 'r') as infile:
        for l in infile:
            a, b, n_kmers, ident, qstrand, qstart, qend, qlen, tstrand, tstart, tend, tlen = l.split()
            id_a = names.get_id(a)
            id_b = names.get_id(b)
            
            # Check if sequences are containd. Do not want contigs that subsets of other contigs.
            # Allow for 10 mismatches at ends
            if int(qend) - int(qstart) < int(qlen) - 10 and int(tend) - int(tstart) < int(tlen) - 10:
                candidate_edges.append((id_a, id_b))
            else:
                # if not quiet:
                #     sys.stderr.write("Dropped overlap between " + a + " and " + b + " because of containment.\n")
                if int(qend) - int(qstart) < int(qlen) - 10: # Want to ignore id_a
                    ignore_vertices[id_a] = True
                else:           # Ignore id_b instead
                    ignore_vertices[id_b] = True
                    
        for (id_a, id_b) in candidate_edges:
            g.add_edge(id_a, id_b)
            
        if g.n_vertices() > 0:
            return g, names, len(ignore_vertices)
        else:
            raise IOError('No vertices in input graph')


def output_components(cl, names, fname, oh):
    '''
    Input: component list, vertex names, input filename, output file
    Output: JSON with list of components, each component being a list of names.
    '''
    oh.write('{\n')
    oh.write('  "graph_file": ' + '"' + fname + '",\n')
    oh.write('  "components": ' + '[\n')
    first = True
    for c in cl:
        if not first:
            oh.write(',\n')
        size = len(c)
        sorted_vertex_list = sorted(list(c))
        s_vertices = translate_component(sorted_vertex_list, names) # List of vertex names
        md5 = hashlib.md5(s_vertices.encode('utf-8')).hexdigest()
        
        oh.write('    { "size"  : ' + str(size) + ',\n      "md5"   : "' + md5 + '",\n      "nodes" : ' + str(s_vertices) + '\n    }')
        first = False
    oh.write('\n  ]\n}\n')
    
def output_components_lst(cl, names, prefix):
    '''
    Input: component list, vertex names, input filename, output file
    Output: a simple list of contig ids, one file per component.
    '''
    i = 0
    for c in cl:
        with open(prefix + str(i) + ".lst", "w") as oh:
            for id in c:
                name = names.get_name(id)
                oh.write(name + "\n")
   

def translate_component(c, names):
    l = []
    for v in c:
        l.append('"' + names.get_name(v) + '"')
    s = "[" + ", ".join(l) + "]"
    return s


def output_degree_distribution(g, oh):
    degree_distr = dict()
    for v in g.vertices():
        d = g.degree(v)
        if d in degree_distr:
            degree_distr[d] += 1
        else:
            degree_distr[d] = 1

    for v, d in degree_distr.items():
        oh.write("{:d}\t{:d}\n".format(v, d))


def filter_components_by_size(cl, minzise, maxsize):
    '''
    If minzise/maxsize set, remove components outside the range.
    '''
    if minzise:
        cl = list(filter(lambda c: len(c) >= minzise, cl))
    if maxsize:
        cl = list(filter(lambda c: len(c) <= maxsize, cl))
    return cl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph_file", help="Input graph in M4 format. See olp_overlap output.")
    parser.add_argument("-g", "--getcomponents", type=str, help="Output components to named file.")
    parser.add_argument("-G", "--degreelimit", type=int, help="Limit components to given degree. High degree vertices can be 'endpoints' of a component, but will not be further traversed.")
    parser.add_argument("-l", "--lst", action="store_true", help="Like option -g, but output in one component per file, one contig accession per line.")
    parser.add_argument("-r", "--randomsubgraphs", type=int, help="Output <int> random subgraphs. Restricted to components of size > 5.")
    parser.add_argument("-m", "--maxdegree", type=int, help="Specify a maximum allowed degree on a vertex. Vertices with higher degree is removed from the graph.")
    parser.add_argument("-cmin", "--mincomponentsize", type=int, help="Forget about components smaller than this.")
    parser.add_argument("-cmax", "--maxcomponentsize", type=int, help="Forget about components larger than this.")
    parser.add_argument("-d", "--degreedistribution", type=str, help="Compute and write the degree distribution (histogram) of the graph to the named output file.")
    parser.add_argument("-s", "--simplestats", action="store_true", help="Just output number of vertices and components.")
    parser.add_argument("-q", "--quiet", action="store_true", help="Do not print progress information.")
    
    args=parser.parse_args()
    
    if not args.quiet:
        sys.stderr.write('Reading graph\n')
    g, names, n_dropped = read_graph(args.graph_file, args.quiet)
    if not args.quiet:
        sys.stderr.write('Ignored ' + str(n_dropped) + ' contained contigs\n')

    if args.maxdegree:
        if not args.quiet:
            sys.stderr.write('Removing high degree nodes\n')
        n = g.remove_high_degree_nodes(args.maxdegree)
        sys.stderr.write("Removed " + str(n) + " high degree nodes\n")
        n = g.remove_singletons()
        sys.stderr.write("Removed " + str(n) + " singletons\n")
        

    if args.degreedistribution:
        with open(args.degreedistribution, "w") as oh:
            output_degree_distribution(g, oh)

    if args.getcomponents or args.lst or args.randomsubgraphs or args.simplestats:
        if not args.quiet:
            sys.stderr.write('Collecting components\n')
        # Get a list of connected components
        if args.degreelimit:
            cl = g.degree_limited_connected_components(args.degreelimit)
        else:
            cl = g.connected_components() 

        cl = list(filter_components_by_size(cl, args.mincomponentsize, args.maxcomponentsize))

        if args.getcomponents:  # JSON format
            with open(args.getcomponents, "w") as oh:
                output_components(cl, names, args.graph_file, oh)

        if args.lst:          # One component per file and plain single id per line
            output_components_lst(cl, names, "component_")


        if args.randomsubgraphs:
            if not args.randomsubgraphs > 0:
                raise Exception('Must ask for a positive number of subgraphs!')
            if len(cl) == 0:
                raise Exception('No components available!')
            for i in range(args.randomsubgraphs):
                r = random.randint(0, len(cl)-1)
                filename='component_r' + str(i) + '.dot'
                g.output_component_dot(g, names, filename, cl[r])
            
        if args.simplestats:
            print("n vertices:   ", g.n_vertices())
            print("n components: ", len(list(cl)))


if __name__ == '__main__':
    main()
