
This is a C++ Wrapper for igraph library. This project’s goal is to
automatically manage the life cycle of igraph objects and to provide
exception mechanism to handle errors.

## List of objects

#### Graph
* Graph — `igraph_t`

#### Vector
* Vector — `igraph_vector_t`
* VectorPtr<> — `igraph_vector_ptr_t`

#### Iterator
* VertexIterator — `igraph_vit_t`
* EdgeIterator — `igraph_eit_t`

#### Selector
* VertexSelector — `igraph_vs_t`
* EdgeSelector — `igraph_es_t`

#### Matrix
* Matrix — `igraph_matrix_t`

#### Sparse matrix
* ~~`igraph_spmatrix_t`~~
* ~~`igraph_sparsemat_t`~~

#### Stack
* ~~`igraph_stack_t`~~

#### Double-ended queues
* ~~`igraph_dqueue_t`~~

#### Maximum and minimum heaps
* ~~`igraph_heap_t`~~

#### String vector
* ~~`igraph_strvector_t`~~

#### Adjacency list
* ~~`igraph_adjlist_t`~~

#### Incident edges
* ~~`igraph_inclist_t`~~
* ~~`igraph_lazy_inclist_t`~~

#### Random numbers
* ~~`igraph_rnd_t`~~

#### Attributes
* ~~`igraph_attribute_table_t`~~

#### Hierarchical random graphs
* ~~`igraph_hrg_t`~~

## Graph generators

### Deterministic graph generators
* Create — `igraph_create`
* Small — `igraph_small`
* Adjacency matrix — `igraph_adacency`
* ~~Weighted adjacency matrix — `igraph_weighted_adjacency`~~
* ~~Adjacency list — `igraph_adjlist`~~
* Star — `igraph_star`
* Lattice — `igraph_lattice`
* Ring — `igraph_ring`
* Tree — `igraph_tree`
* Full — `igraph_full`
* Full citation — `igraph_full_citation`
* Famous — `igraph_famous`
* LCF — `igraph_lcf` & `igraph_lcf_vector`
* Atlas — `igraph_atlas`
* de Bruijn — `igraph_de_bruijn`
* Kautz — `igraph_kautz`
* Extended chordal ring — `igraph_extended_chordal_ring`
* Connect neighborhood — `igraph_connect_neighborhood`

### Randomized graph generators
* ~~Geometric random  — `igraph_grg_game`~~
* ~~Barabási–Albert — `igraph_barabasi_game`~~
* Erdős–Rényi — `igraph_erdos_renyi_game`
* ~~Watts–Strogatz — `igraph_watts_strogatz_game`~~
* ~~Rewire edges — `igraph_rewire_edges`~~
* ~~Degree sequence — `igraph_degree_sequence_game`~~
* ~~k-regular  — `igraph_k_regular_game`~~
* ~~Static fitness — `igraph_static_fitness_game`~~
* ~~Static power-law — `igraph_static_power_law_game`~~
* ~~Forest fire — `igraph_forest_fire_game`~~
* ~~Rewire  — `igraph_rewire`~~
* ~~Growing random — `igraph_growing_random_game`~~
* ~~Callaway traits — `igraph_callaway_traits_game`~~
* ~~Establishment — `igraph_establishment_game`~~
* ~~Preference — `igraph_preference_game`~~
* ~~Asymmetric preference — `igraph_asymmetric_preference_game`~~
* ~~Recent degree — `igraph_recent_degree_game`~~
* ~~Barabási aging — `igraph_barabasi_aging_game`~~
* ~~Recent degree aging — `igraph_recent_degree_aging_game`~~
* ~~Cited type — `igraph_cited_type_game`~~
* ~~Citing cited type — `igraph_citing_cited_type`~~
* ~~Stochastic block model — `igraph_sbm_game`~~

## Games on graphs

### Microscopic update rules
* ~~`igraph_deterministic_optimal_imitation`~~
* ~~`igraph_moran_process`~~
* ~~`igraph_roulette_wheel_imitation`~~
* ~~`stochastic_imitation`~~

## Library functions
* Structural properties — *Partial*
* ~~Neighborhood~~
* ~~Graph components~~
* ~~Degree sequences~~
* ~~Centrality measures~~
* ~~Estimating centality measures~~
* ~~Centralization~~
* ~~Similarity measures~~
* ~~Spanning trees~~
* ~~Transitivity and clustering coefficient~~
* ~~Directedness conversion~~
* ~~Spectral properties~~
* ~~Non-simple graphs: multiple and loop edges~~
* ~~Mixing patterns~~
* ~~K-Cores~~
* ~~Topological sorting & directed acyclic graphs~~
* ~~Maximum cardinality search, graph decomposition, chordal graphs~~
* ~~Matchings~~
* ~~Line graphs~~
* ~~Unfolding a graph into a tree~~
* ~~Other operations~~
* ~~Graph visitors, BFS & DFS~~
* ~~Cliques & Independent Vertex Sets~~
* ~~Graph Isomorphism~~
* ~~Motifs, dyad census & triad census~~
* ~~Graph layouts~~
* ~~Maximum flows, minimum cuts & releated measures~~
* ~~Reading & writing graphs from and to files~~
* ~~Graph layouts~~
* ~~Vertex separatos~~
* ~~Community structure~~
* ~~Graphlets~~
* ~~Spectral coase graining~~
* ~~Graph operators~~
* ~~BLAS~~
* ~~Bipartite~~
* ~~Multi-thread support~~
* ~~Not graph-related functions~~













