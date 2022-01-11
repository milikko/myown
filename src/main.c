#include "main.h"

/*
  indicates whether the algorithms works with 1 or 0 indexed data.
  If e.g. 1-indexed, the offset is 1 such that the backend (so how the data
  is stored) is 0-indexed. But ALL print out would still be 1-indexed as desired.
  The offset is measured in the input data, by finding the lowest value for the vertex id
 */
int offset = 0;

/**
 * help - a help utility
 * writing out a help utility to the user of the program
 * to be called with flag --help when executing the program
 */
void help () {
	printf ("To run, required arguments are as follows (separate each argument with a whitespace):\n\n");
	printf ("./main <algorithm> <inputfile> <outputfile> <k integer> <u integer> <v integer> <query times>\n");
	printf ("\nPossible input for each flag:\n\n");
	printf ("<algorithm>: dj, djopt, tz, bdj\n");
	printf ("<inputfile>: tests/USANY.txt (needs to be of DIMACS ssp file format) \n");
	printf ("<outputfile>: output.csv (will be generated automatically) \n");
	printf ("<k integer>: Any integer that satisfies k >= 1 \n");
	printf ("<u integer>: Any integer that represents a node from <inputfile> \n");
	printf ("<v integer>: Any integer that represents a node from <inputfile> \n");
	printf ("<query times>: Any positive integer, indicating how many times the Thorup-Zwick query alg. shall be executed\n\n");
}

/**
 * run_bdj - wrapper function for running bidirectional dijkstra from u to v
 * @graph: graph with vertices and edges
 * @u: source vertex u
 * @v: target vertex v
 * @n: number of vertices in graph
 * @m: number of edges in graph
 * Calls Bidirectional Dijkstra's algorithm, and measures the spent RAM and CPU time
 */
struct ssp_res *run_bdj (struct graph *graph, int u, int v, int n, int m)
{
	struct ssp_res *bdj;
	double cpu_time_spent = 0.0, avg_ins = 0.0, avg_dec = 0.0, avg_ext = 0.0;
	for (int i = 0; i < DIJKSTRA_TIMES; i++) {
		clock_t begin = clock();
		bdj = bidirectional_dijkstra (graph, u-offset, v-offset);
		clock_t end = clock();
		if (bdj == NULL) {
			perror ("No path (u, v) in djkstra_opt_alg could be found\n");
			exit (-1);
		}
		cpu_time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
		avg_ext += bdj->avg_extract_min_time /
			bdj->num_extract_min / MS;
		avg_dec += bdj->avg_decrease_key_time /
			bdj->num_decrease_key / MS;
		avg_ins += bdj->avg_min_heap_insert_time = bdj->avg_min_heap_insert_time /
			bdj->visited_nodes / MS;
	}
	bdj->memory_consump = get_vm_peak ();
	bdj->avg_extract_min_time = (avg_ext / (double)DIJKSTRA_TIMES);
	bdj->avg_decrease_key_time = (avg_dec / (double)DIJKSTRA_TIMES);
	bdj->avg_min_heap_insert_time = (avg_ins / (double)DIJKSTRA_TIMES);
	bdj->dist_time = (cpu_time_spent / DIJKSTRA_TIMES);
	bdj->visited_nodes_ratio = (double)bdj->visited_nodes/(double)(2*n);
	bdj->visited_edges_ratio = (double)bdj->visited_edges/(double)(2*m);

	bdj->query_times = DIJKSTRA_TIMES;

	printf ("\nResult of Bidirectional Dijkstra (%d, %d) = %d\n", u, v, bdj->dist);
	printf ("vertices n=%d, edges m=%d\n", n, m);
	printf ("Visiting ratio of vertices = %f, edges = %f\n",
			bdj->visited_nodes_ratio,
			bdj->visited_edges_ratio);
	printf ("%d extract-min operations. Avg time pr. operation: %.10f ms\n",
			bdj->num_extract_min, bdj->avg_extract_min_time);
	printf ("%d decrease-key operations. Avg time pr. operation: %.10f ms\n",
			bdj->num_decrease_key, bdj->avg_decrease_key_time);
	printf ("%d min-heap-insert operations. Avg time pr. operation: %.10f ms\n",
			bdj->visited_nodes, bdj->avg_min_heap_insert_time);
	printf ("Time spent on running Bidirectional Dijkstra (%d, %d) = %f sec\n", u, v, bdj->dist_time);
	printf ("Algorithm is executed %d times\n", bdj->query_times);
	printf ("Memory usage of Bidirectional Dijkstra = %d KB\n", bdj->memory_consump);

	return bdj;
}

/**
 * run_dijkstra - wrapper function for running dijkstra from source vertex
 * @graph: graph with vertices and edges
 * @u: source vertex u
 * @v: target vertex v
 * @n: number of vertices in graph
 * @m: number of edges in graph
 * Calls Dijkstra's algorithm, and measures the spent RAM and CPU time
 */
struct ssp_res *run_dijkstra (struct graph *graph, int u, int v, int n, int m)
{
	struct ssp_res *dijkstra;
	double cpu_time_spent = 0.0, avg_dec = 0.0, avg_ext = 0.0;
	for (int i = 0; i < DIJKSTRA_TIMES; i++) {
		clock_t begin = clock();
		dijkstra = dijkstra_alg (graph, u-offset);
		clock_t end = clock();
		if (dijkstra == NULL) {
			perror ("No path (u, v) in djkstra_opt_alg could be found\n");
			exit (-1);
		}
		cpu_time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
		avg_ext += dijkstra->avg_extract_min_time /
			dijkstra->visited_nodes / MS;
		avg_dec += dijkstra->avg_decrease_key_time /
			dijkstra->num_decrease_key / MS;
	}
	dijkstra->memory_consump = get_vm_peak ();
	dijkstra->avg_extract_min_time = (avg_ext / (double)DIJKSTRA_TIMES);
	dijkstra->avg_decrease_key_time = (avg_dec / (double)DIJKSTRA_TIMES);
	dijkstra->dist_time = (cpu_time_spent / DIJKSTRA_TIMES);
	dijkstra->visited_nodes_ratio = (double)dijkstra->visited_nodes/(double)n;
	dijkstra->visited_edges_ratio = (double)dijkstra->visited_edges/(double)m;

	dijkstra->dist = dijkstra->S_f[v-offset].sp_est;

	dijkstra->query_times = DIJKSTRA_TIMES;

	printf ("\nResult of Dijkstra's algorithm (%d, %d) = %d\n", u, v,
			dijkstra->dist);
	printf ("vertices n=%d, edges m=%d\n", n, m);
	printf ("Visiting ratio of vertices = %f, edges = %f\n",
			dijkstra->visited_nodes_ratio,
			dijkstra->visited_edges_ratio);
	printf ("%d extract-min operations. Avg time pr. operation: %.10f ms\n",
			dijkstra->num_extract_min, dijkstra->avg_extract_min_time);
	printf ("%d decrease-key operations. Avg time pr. operation: %.10f ms\n",
			dijkstra->num_decrease_key, dijkstra->avg_decrease_key_time);
	printf ("%d min-heap-insert operations. Avg time pr. operation: %.1f ms\n",
			0, 0.0);
	printf ("Time spent on running Dijkstra's algorithm (%d, %d) = %f sec\n",
			u, v, dijkstra->dist_time);
	printf ("Algorithm is executed %d times\n", dijkstra->query_times);
	printf ("Memory usage of Dijkstra's algorithm = %d KB\n", dijkstra->memory_consump);

	return dijkstra;
}

/**
 * run_dijkstra - wrapper function for running dijkstra from source vertex
 * @graph: graph with vertices and edges
 * @u: source vertex u
 * @v: target vertex v
 * @n: number of vertices in graph
 * @m: number of edges in graph
 * Calls Dijkstra's algorithm, and measures the spent RAM and CPU time
 */
struct ssp_res *run_opt_dijkstra (struct graph *graph, int u, int v, int n, int m)
{
	struct ssp_res *dijkstra;
	double cpu_time_spent = 0.0, avg_ins = 0.0, avg_dec = 0.0, avg_ext = 0.0;
	for (int i = 0; i < DIJKSTRA_TIMES; i++) {
		clock_t begin = clock();
		dijkstra = dijkstra_opt_alg (graph, u-offset, v-offset);
		clock_t end = clock();
		if (dijkstra == NULL) {
			perror ("No path (u, v) in dijkstra_opt_alg could be found\n");
			exit (-1);
		}
		cpu_time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
		avg_ins += dijkstra->avg_min_heap_insert_time /
			dijkstra->visited_nodes / MS;
		avg_ext += dijkstra->avg_extract_min_time /
			dijkstra->num_extract_min / MS;
		avg_dec += dijkstra->avg_decrease_key_time /
			dijkstra->num_decrease_key / MS;
	}
	dijkstra->memory_consump = get_vm_peak ();
	dijkstra->avg_min_heap_insert_time = (avg_ins / (double)DIJKSTRA_TIMES);
	dijkstra->avg_extract_min_time = (avg_ext / (double)DIJKSTRA_TIMES);
	dijkstra->avg_decrease_key_time = (avg_dec / (double)DIJKSTRA_TIMES);
	dijkstra->dist_time = (cpu_time_spent / DIJKSTRA_TIMES);
	dijkstra->visited_nodes_ratio = (double)dijkstra->visited_nodes/(double)n;
	dijkstra->visited_edges_ratio = (double)dijkstra->visited_edges/(double)m;

	dijkstra->dist = dijkstra->S_f[v-offset].sp_est;

	dijkstra->query_times = DIJKSTRA_TIMES;

	printf ("\nResult of optimised Dijkstra's algorithm (%d, %d) = %d\n", u, v,
			dijkstra->dist);
	printf ("vertices n=%d, edges m=%d\n", n, m);
	printf ("Visiting ratio of vertices = %f, edges = %f\n",
			dijkstra->visited_nodes_ratio, dijkstra->visited_edges_ratio);
	printf ("%d extract-min operations. Avg time pr. operation: %.10f ms\n",
			dijkstra->num_extract_min, dijkstra->avg_extract_min_time);
	printf ("%d decrease-key operations. Avg time pr. operation: %.10f ms\n",
			dijkstra->num_decrease_key, dijkstra->avg_decrease_key_time);
	printf ("%d min-heap-insert operations. Avg time pr. operation: %.10f ms\n",
			dijkstra->visited_nodes, dijkstra->avg_min_heap_insert_time);
	printf ("Time spent on running optimised Dijkstra (%d, %d) = %f sec\n",
			u, v, dijkstra->dist_time);
	printf ("Algorithm is executed %d times\n", dijkstra->query_times);
	printf ("Memory usage of optimised Dijkstra = %d KB\n", dijkstra->memory_consump);

	return dijkstra;
}

#include "thorupzwick.h"

void pp_aseqs (struct aseq **A, int k)
{
	printf("Pretty print all A_i sequences\n");
	for (int i = 0; i < k; i++) {
		if (A[i]->setlength > 0) {
			for (int j = 0; j < A[i]->setlength; j++)
				printf ("A_{%d,%d} v_id:%d spest:%d\n", i, j,
						A[i]->nodes[j].v_id+offset, A[i]->nodes[j].sp_est);
		}
	}
	return;
}

void pp_clusters (struct clusterlist **C, int i)
{
	printf("Pretty print all clusters\n");
	for (int a = 0; a < C[i]->num_clusters; a++)
		for (int b = 0; b < C[i]->clusters[a].num_nodes; b++)
			printf ("cluster nodes for C(w:%d)=v:%d, SP:%d\n",
					C[i]->clusters[a].w->v_id+offset,
					C[i]->clusters[a].nodes[b].v_id+offset,
					C[i]->clusters[a].nodes[b].sp_est);

	return;
}

void pp_bunches (struct bunchlist *bunchlist, int k)
{
	printf("Pretty print all bunches\n");
	for (int i = 0; i < bunchlist->num_bunches; i++) {
		struct node *s;
		// iterate over hash elements of B(v)
		for (s = bunchlist->bunches[i].nodes; s != NULL; s = s->hh.next) {
			printf("B(v:%d)=w:%d, SP:%d\n",
				   bunchlist->bunches[i].v->v_id+offset,
				   s->v_id+offset, s->sp_est);
		}
		for (int j = 0; j < k; j++) {
			printf("p_%d(%d)=%d, SP:%d\n",
				   j, bunchlist->bunches[i].v->v_id+offset,
				   bunchlist->bunches[i].piv[j].v_id+offset,
				   bunchlist->bunches[i].piv[j].sp_est);
		}
	}
	return;
}

void pp_pivots (struct bunchlist *bunchlist, struct node *nodes,
				unsigned int n, int i)
{
	printf("Pretty print all pivot elements (witnesses) \n");
	for (unsigned int j = 0; j < n-1; j++) {
		printf ("p_%d(%d) = %d, dist:%d\n", i, nodes[j].v_id+offset,
				bunchlist->bunches[nodes[j].v_id].piv[i].v_id+offset,
				bunchlist->bunches[nodes[j].v_id].piv[i].sp_est);
	}
}

/*
 ************************************************************************
 * Begin of the actual functions for thorup-zwick
 ************************************************************************
*/

/**
 * add_s_node_to_graph - adding the super vertex s to the graph
 * @graph: graph
 * @ai: all sequences
 * The added node s will always have the vertex id of the current
 * number of vertices in the graph
 */
void add_s_node_to_graph (struct graph *graph, struct aseq *ai)
{
	for (int i = 0; i < ai->setlength; i++) {
		add_edges (graph, graph->V, ai->nodes[i].v_id, 0);
	}
	return;
}

/**
 * initialise_single_source_tz - initialising heap for Dijkstra's alg.
 * @graph: for the graph pointers
 * Running time: O(n)
 * This copies pointers of heap nodes to graph, to keep bidirectional
 * pointers. Also it initialises a heap but with no single source actually
 * - this is set manually outside the function
 */
struct heap *initialise_single_source_tz (unsigned int n)
{
	struct heap *heap = malloc (sizeof (struct heap));
	// One extra slot for node s (later used for dijkstra)
	heap->nodes = malloc ((n+1) * sizeof(struct node*));

	if (heap == NULL || heap->nodes == NULL) {
		perror ("Pointer error of heap\n");
		exit (-1);
	}

	int val = (int) INFINITY;

	for (unsigned int i = 0; i < n; i++) {
		heap->nodes[i] = add_node (i, val, i);
	}

	heap->heap_size = n;

	return heap;
}

/**
 * dijkstra_alg_tz - Dijkstra's algorithm (modified for Thorup-Zwick)
 * @graph: the graph G = (V, E)
 * @Q: the heap to work with.
 * Running time: O((m + n) lg n)
 * Executing dijkstra's algorithm but without a specified single source as
 * initialise single source procedure is required to be called on before hand
 * The algorithm is implemented with a priorityqueue (min-heap) as a data structure.
 * pop takes lg n, V = E = m = n, loop m times
 * Is used to for computing pivots
 */
struct node *dijkstra_alg_tz (struct graph *graph, struct heap *Q)
{
	struct node *S = malloc (Q->heap_size * sizeof (struct node));

	if (S == NULL) {
		perror ("Pointer error of S array\n");
		exit (-1);
	}

	while (Q->heap_size != 0) {
		struct node *u = extract_min (Q);
		memcpy (&S[u->v_id], u, sizeof (struct node));
		for (struct adjlistnode *s = graph->adjlists[u->v_id].head;
			 s != NULL; s = s->next) {
			struct node *v = graph->adjlists[s->v_id].nd;
			if ((v != NULL) && (v->sp_est > u->sp_est + s->weight)) {
				int sp_est = u->sp_est + s->weight;
				decrease_key (Q, v, u, sp_est);
			}
		}
		// avoiding the u node can be relaxed again
		graph->adjlists[u->v_id].nd = NULL;
	}

	free_heap (Q);
	return S;
}

/**
 * find_pivot - finding all pivot elements for v in A_i
 * @aiplusone_pivot_arr: A_i+1 set
 * @n: number of nodes in graph, |V|=n
 * Running time: O(m)
 * Not looping through all v in V and A sets, otherwise running time could be O(n^2)
 * If v.pi = s, it is a pivot node
 * If a node has been visited by a traversing one, save the visitor, such that
 * in a later iteration we check if the respective node has been visited.
 * If so, we simply retrieve the visitor node's pivot node
 * (*) property, if d(A_i, v) == d(A_i+1,v), copy! Means we get p_i(v) = p_i+1(v) ∈ B(v),
 * otherwise, δ(p_i(v), v) = δ(A_i, v) < δ(A_i+1, v). Thus distance will be same,
 * but the node would change (to A_i+1) if (*) fails
 */
struct node *find_pivot (struct node *aiplusone_pivot_arr,
						 struct node *nodes, unsigned int n)
{
	// no need to allocate for s node, therefore n-1
	int *visited_nodes = malloc ((n-1) * sizeof (int));
	for (unsigned int i = 0; i < (n-1); i++)
		visited_nodes[i] = -1;

	struct node *pivot_arr = malloc ((n-1) * sizeof (struct node));
	memset (pivot_arr, 0, ((n-1) * sizeof (int)));

	if (pivot_arr == NULL) {
		perror ("Failed to allocate pivot array\n");
		exit (-1);
	}

	// skipping node s in nodes, assuming s is the last node
	for (unsigned int i = 0; i < n-1; i++) {
		struct node *piv = NULL;
		struct node *v = nodes[i].pi;

		if (visited_nodes[nodes[i].v_id] != -1) {
			memcpy (&pivot_arr[nodes[i].v_id], &pivot_arr[visited_nodes[nodes[i].v_id]],
					sizeof(pivot_arr[visited_nodes[nodes[i].v_id]]));
			pivot_arr[nodes[i].v_id].sp_est = nodes[nodes[i].v_id].sp_est;
		} else if (v->v_id == nodes[n-1].v_id) {
			memcpy (&pivot_arr[nodes[i].v_id], &nodes[i], sizeof(nodes[i]));
			pivot_arr[nodes[i].v_id].sp_est = 0;
		} else {
			while (v) {
				piv = v;
				// avoid reaching the s node
				if (v->pi != NULL && v->pi->v_id == nodes[n-1].v_id) {
					break;
				}
				visited_nodes[v->v_id] = nodes[i].v_id;
				v = v->pi;
			}
			// let piv be w, thus we set p_i(v) <- w
			memcpy (&pivot_arr[nodes[i].v_id], piv, sizeof(&piv));
			pivot_arr[nodes[i].v_id].sp_est = nodes[nodes[i].v_id].sp_est;
		}
		// (*) property
		if (aiplusone_pivot_arr != 0 && pivot_arr[nodes[i].v_id].sp_est
			== aiplusone_pivot_arr[nodes[i].v_id].sp_est) {
			memcpy (&pivot_arr[nodes[i].v_id], &aiplusone_pivot_arr[nodes[i].v_id],
					sizeof (aiplusone_pivot_arr[nodes[i].v_id]));
		}
	}

	FREE (visited_nodes);

	return pivot_arr;
}

/**
 * dijkstra_cluster_tz - modified version of dijkstra for constructing clusters
 * @graph: graph
 * @w: all sequences
 * @pivot_nodes: all pivot nodes
 * @k: k integer
 * @i: the current i index
 * Running time: O(m)
 * Cluster C(w) is composed of all the vertices that are closer to w than any (i+1)-center
 * Note that δ(A_{i+1}, v) was computed in the previous iteration so the test takes constant time.
 * Note that if v not in C(w), then we never assign a finite distance, as by definition,
 * v ∈ C(w) if and only if δ(w, v) < δ(A_i+1, v). This modified version of Dijkstra starts with w in the heap,
 * consequently adding all vertices being relaxed (or relaxed if already added to heap).
 * We relax the edge (u,v) ONLY if u.d + l(u,v) < d(A_{i+1}, v), where d(A_{i+1}, v) = d(p_{i+1}(v), v)
 * Also, if i==k-1, no pivot element to compare with, thus we just relax edge / insert node
 */
struct cluster *dijkstra_cluster_tz (struct graph *graph, struct node *w,
									 struct node *pivot_nodes, int k, int i)
{
	int idx = 0;
	int num_nodes = 0;
	int *in_heap = calloc (graph->V, sizeof(int));
	int *extracted = calloc (graph->V, sizeof(int));
	struct cluster *cluster = malloc (sizeof (struct cluster));
	struct heap *Q = malloc (sizeof (struct heap));
	Q->nodes = malloc (graph->V * sizeof(struct node*));
	struct node *S = malloc (graph->V * sizeof (struct node));

	if (Q == NULL || S == NULL || in_heap == NULL
		|| extracted == NULL || Q->nodes == NULL || cluster == NULL) {
		perror ("Failed to allocate memory for in_heap, extracted, S, Q, Q->nodes or cluster\n");
		exit (-1);
	}

	Q->heap_size = 0;
	// w.d = 0
	min_heap_insert (Q, w->v_id, 0, graph);
	num_nodes += 1;

	while (Q->heap_size != 0) {
		struct node *u = extract_min (Q);
		// tmp_S[idx] gets assigned value from u, copying to its address
		memcpy (&S[idx], u, sizeof(struct node));
		// avoiding node u can be relaxed in next iterations,
		// as u is no longer on heap
		extracted[u->v_id] = 1;
		in_heap[u->v_id] = 0;
		for (struct adjlistnode *s = graph->adjlists[S[idx].v_id].head;
			 s != NULL; s = s->next) {
			struct node *v = graph->adjlists[s->v_id].nd;
			if ((i == k-1 && extracted[s->v_id] == 0) ||
				((pivot_nodes != 0) &&
					 (u->sp_est + s->weight) < pivot_nodes[s->v_id].sp_est
					 && extracted[s->v_id] == 0)) {
				int sp_est = u->sp_est + s->weight;
				if (in_heap[s->v_id] == 1) {
					if ((v->sp_est > sp_est)) {
						decrease_key (Q, v, u, sp_est);
					}
				} else {
					min_heap_insert (Q, s->v_id,
									 sp_est, graph);
					graph->adjlists[s->v_id].nd->pi = u;
					in_heap[s->v_id] = 1;
					num_nodes += 1;
				}
			}
		}
		idx += 1;
	}

	cluster->nodes = malloc (num_nodes * sizeof (struct node));
	if (cluster->nodes == NULL) {
		perror ("Alloc of cluster->nodes failed\n");
		exit (-1);
	}
	for (int i = 0; i < num_nodes; i++)
		memcpy (&cluster->nodes[i], &S[i], sizeof (struct node));

	cluster->w = w;
	cluster->num_nodes = num_nodes;

	FREE (extracted);
	FREE (in_heap);
	FREE (S);
	free_heap (Q);

	return cluster;
}

/**
 * construct_clusters - constructing clusters for all w \in A_i - A_i+1
 * @graph: adjacency list rep. of the graph
 * @A: all sequences
 * @pivot_nodes: all pivot nodes
 * @i: the current i index
 * @k: k integer
 * Running time: O(n)
 * From here we construct the clusters. If i=k-1 the clusters of all w in A_i contain
 * all vertices. Otherwise, we check if w is in A_i - A_i+1. In both cases we compute
 * a modified version of dijkstra from with w as source.
 */
struct clusterlist *construct_clusters (struct graph *graph, struct aseq **A,
										struct node *pivot_nodes, int i, int k)
{
	int num_clusters = 0;
	struct cluster *tmp_clusters[A[i]->setlength];
	struct clusterlist *C = malloc (sizeof (struct clusterlist));

	// For every w in A_i - A_i+1
	for (int w = 0; w < A[i]->setlength; w++) {
		// for every w ∈ A_k−1 we have C(w) = V
		// if not i==k-1, check if w in A_i - A_i+1
		if (i < (k-1) && A[i+1]->added[A[i]->nodes[w].v_id]) {
			continue;
		}
		// if indeed w in A_i - A_i+1, go on with dijkstra to create a cluster C(w)
		// For i==k-1, create cluster C(w) = V,
		// computing the distances from w to all other v in V
		struct cluster *cw = dijkstra_cluster_tz (graph, &A[i]->nodes[w], pivot_nodes, k, i);
		tmp_clusters[num_clusters] = malloc (sizeof (struct cluster));
		memcpy (tmp_clusters[num_clusters], cw, sizeof (struct cluster));
		num_clusters += 1;
	}

	C->clusters = malloc (num_clusters * sizeof (struct cluster));
	if (C->clusters == NULL) {
		perror ("Alloc of C->clusters failed\n");
		exit (-1);
	}
	for (int i = 0; i < num_clusters; i++) {
		memcpy (&C->clusters[i], tmp_clusters[i], sizeof (struct cluster));
		FREE (tmp_clusters[i]);
	}

	C->num_clusters = num_clusters;

	return C;
}

/**
 * compute_bunch_size - computing the size of the bunches
 * @bunchlist: adjacency list rep. of the graph
 * Running time: O(n)
 * Computes the size of all bunches by computing the hash tables overhead.
 */
size_t compute_bunch_size (struct bunchlist *bunchlist)
{
	size_t ds_size = 0;
	for (int i = 0; i < bunchlist->num_bunches; i++) {
		struct node *s = bunchlist->bunches[i].nodes;
		ds_size += HASH_OVERHEAD (hh, s);
	}
	return ds_size;
}

/**
 * construct_bunches - constructing bunches for all v in V
 * @C: all clusters
 * @k: k integer
 * @A: all sequences
 * @bunchlist: a memory initialised struct of bunches
 * @pivot_nodes: all pivot nodes
 * Running time: Runs in linear time of the cluster size
 * such that O(\sum_{w ∈ V} |C(w)|) = O(\sum_{v ∈ V} |B(v)|)
 * Generating bunches from the clusters. For each v in V, we construct a bunch,
 * containing all w of all C(w) v belongs to. Note, we get A_{k-1} \subseteq B(v)
 * for every v in V, meaning all nodes from A_{k-1} belong to every bunch
 */
void construct_bunches (struct clusterlist **C, int k,
						struct bunchlist *bunchlist)
{
	for (int i = 0; i < k; i++) {
		for (int c = 0; c < C[i]->num_clusters; c++) {
			struct cluster *cluster = &C[i]->clusters[c];
			for (int v = 0; v < cluster->num_nodes; v++) {
				struct node *s;
				int tmp = cluster->nodes[v].v_id;
				HASH_FIND_INT (bunchlist->bunches[tmp].nodes, &cluster->w->v_id, s);
				if (s == NULL) {
					cluster->w->sp_est = cluster->nodes[v].sp_est;
					s = malloc (sizeof (struct node));
					memcpy (s, cluster->w, sizeof(&cluster->w));
					// adding node to the hash
					HASH_ADD_INT (bunchlist->bunches[tmp].nodes, v_id, s);
					bunchlist->bunches[tmp].num_nodes += 1;
				}
			}
		}
	}

	bunchlist->bunch_size += compute_bunch_size (bunchlist);

	return;
}

/**
 * create_aseqs - constructing sets of vertices for A_0, ... A_k-1
 * @A: the A sequences to fill in
 * @k: k integer
 * @graph: adjacency list rep. of graph
 * @nodes: all n nodes and their attributes
 * Running time: O(n)
 * Constructs non increasing sequence of sets A_0 ⊇ A_1 ⊇ ... ⊇ A_k−1, A_k = Ø
 * Begins with A_0 = V and iterates then 1 <= i < k. Each node in A_i-1 is tested
 * with n^-1/k probability whether it is obtained in A_i. Vertices in an A_i set
 * are also called i-centers
 */
bool create_aseqs (struct aseq **A, int k, struct graph *graph, struct node *nodes)
{
	int *tmp = malloc (graph->V * sizeof(int));
	// A_k <- Ø
	A[k] = NULL;

	// A_0 <- V. Running time O(n)
	A[0] = malloc (sizeof (struct aseq));
	A[0]->nodes = malloc (graph->V * sizeof(struct node));
	A[0]->added = calloc (graph->V, sizeof(int));
	for (unsigned int i = 0; i < graph->V; i++) {
		memcpy (&A[0]->nodes[i], &nodes[i], sizeof(struct node));
		A[0]->added[i] = 1;
	}
	A[0]->setlength = graph->V;

	for (int i = 1; i <= k-1; i++) {
		A[i] = malloc (sizeof (struct aseq));
		A[i]->added = calloc (graph->V, sizeof(int));
		A[i]->setlength = 0;
		for (int j = 0; j < A[i-1]->setlength; j++) {
			// Generates 0.0 - 1.0
			double rnd = (double)randombytes_uniform(RAND_MAX)/RAND_MAX;
			// Check random number is <= n^{-1/k}
			if (rnd <= pow (graph->V, -1.0/(double)k)) {
				int v_id = A[i-1]->nodes[j].v_id;
				tmp[A[i]->setlength] = v_id;
				A[i]->added[v_id] = 1;
				A[i]->setlength += 1;
			}
		}

		if (A[i]->setlength == 0) {
			FREE (tmp);
			for (int j = 0; j < i; j++) {
				FREE (A[j]->nodes);
				FREE (A[j]->added);
				FREE (A[j]);
			}
			return false;
		}

		A[i]->nodes = malloc (A[i]->setlength * sizeof(struct node));
		for (int l = 0; l < A[i]->setlength; l++) {
			memcpy (&A[i]->nodes[l], &nodes[tmp[l]], sizeof (struct node));
		}
	}

	FREE (tmp);

	return true;
}

struct prepro *prepro2 (struct graph *graph, int k)
{
	unsigned int n = graph->V;
	int val = (int) INFINITY;
	bool empty;
	struct node *nodes = malloc (n * sizeof (struct node));
	struct aseq *A[k+1];
	// k+1, for the kth set where d(A_k, v) = infinity for all v in V
	struct node *dist[k+1];
	struct node *pivot_nodes[k+1];
	struct heap *heap;
	struct clusterlist *C[k];
	struct bunchlist *bunchlist = malloc (sizeof (struct bunchlist));
	struct prepro *prepro = malloc (sizeof (struct prepro));

	bunchlist->num_bunches = n;
	bunchlist->bunches = malloc (bunchlist->num_bunches * sizeof(struct bunch));
	bunchlist->bunch_size = 0;

	// saving all v in V in array nodes and Preparing for constructing the bunches
	for (unsigned int i = 0; i < n; i++) {
		memcpy (&nodes[i], add_node (i, val, i), sizeof (struct node));
		bunchlist->bunches[i].v = malloc (sizeof (struct node));
		memcpy (bunchlist->bunches[i].v, &nodes[i], sizeof(struct node));
		bunchlist->bunches[i].piv = malloc ((k+1) * sizeof (struct node));
		memset (&bunchlist->bunches[i].piv[k], 0, sizeof (struct node));
		// Needed for hashing later, thus NULL
		bunchlist->bunches[i].nodes = NULL;
		bunchlist->bunches[i].num_nodes = 0;
		bunchlist->bunch_size += (sizeof (struct node)) + (k * sizeof (struct node));
	}

	// Creating all A_i sequences
	empty = create_aseqs (A, k, graph, nodes);

	// also, since A_{k-1} = Ø, rerun!
	if (!empty) {
		FREE (nodes);
		FREE (bunchlist->bunches->v);
		FREE (bunchlist->bunches->piv);
		FREE (bunchlist->bunches);
		FREE (bunchlist);
		prepro->nodes = NULL;
		prepro->bunchlist = NULL;
		prepro->success = false;
		prepro->k = -1;
		return prepro;
	}

	// d(A_k, v) = infinity, thus NULL
	dist[k] = NULL;
	// p_k(v) undefined as A_k = Ø
	pivot_nodes[k] = NULL;

	// k iterations
	for (int i = k-1; i >= 0; i--) {
		// initialising the heap for current i
		heap = initialise_single_source_tz (graph->V);
		// copy of graph to work with for current i
		struct graph *write_graph = copy_graph_struct (graph, heap);
		// adding source vertex s to G, weight 0 to all other vertices in A_i
		add_s_node_to_graph (write_graph, A[i]);
		// adding s node to heap as well to allow computing sp
		min_heap_insert (heap, write_graph->V, 0, write_graph);
		write_graph->V += 1;
		n = write_graph->V;

		// running dijkstra once for each i
		dist[i] = dijkstra_alg_tz (write_graph, heap);
		// compute d(A_i, v), finding the pivot elements
		pivot_nodes[i] = find_pivot (pivot_nodes[i+1], dist[i], n);
		// Copy the pivot nodes to the piv of bunches
		for (unsigned int j = 0; j < (n-1); j++) {
			bunchlist->bunches[j].piv[i] = pivot_nodes[i][j];
		}

		C[i] = construct_clusters (graph, A, pivot_nodes[i+1], i, k);

		free_graph (write_graph);
		FREE (dist[i]);
	}

	construct_bunches (C, k, bunchlist);

	prepro->nodes = nodes;
	prepro->bunchlist = bunchlist;
	prepro->success = true;
	prepro->k = k;
	prepro->pivot_nodes = pivot_nodes;

	for (int i = 0; i < k; i++) {
		//FREE (pivot_nodes[i]);
		FREE (A[i]->added);
		FREE (A[i]->nodes);
		FREE (A[i]);
		FREE (C[i]->clusters);
		FREE (C[i]);
	}

	// if no A* uncomment
	/* free_graph (graph); */

	return prepro;
}

/**
 * run_tz - wrapper function for running Thorup-Zwick
 * @graph: graph with vertices and edges
 * @k: k integer
 * @u: source vertex u
 * @v: target vertex v
 * Calls Thorup-Zwick algorithm, and measures the spent RAM and CPU time
 * First it calls prepro, the preprocessing algorithm and then calls dist,
 * the query algorithm
 */
struct tz_res *run_tz2 (struct graph *graph, int k, int u, int v, int n, int m, int query_times)
{
	struct tz_res *tz = malloc (sizeof (struct tz_res));
	struct prepro *pp = malloc (sizeof (struct prepro));
	clock_t begin, end;

	tz->dist = 0, tz->dist_time = 0.0;
	tz->query_times = query_times;

	begin = clock();
	pp->success = false;
	while (!pp->success) {
	  pp = prepro2 (graph, k);
	}
	end = clock();
	tz->prepro_time = (double)(end - begin) / CLOCKS_PER_SEC;
	tz->prepro_memory_consump = get_vm_peak();

	for (int i = 0; i < tz->query_times; i++) {
		begin = clock();
		tz->dist += dist (&pp->nodes[u-offset], &pp->nodes[v-offset], pp->bunchlist);
		end = clock();
		tz->dist_time += (double)(end - begin) / CLOCKS_PER_SEC;
	}

	tz->dist = tz->dist / tz->query_times;
	tz->dist_time = tz->dist_time / tz->query_times;
	tz->dist_memory_consump += pp->bunchlist->bunch_size / 1000;
	tz->k = k;

	printf ("Time spent on prepro k=%d Thorup-Zwick: %f\n", tz->k, tz->prepro_time);
	printf ("vertices n=%d, edges m=%d\n", n, m);
	printf ("Memory usage of prepro = %d KB\n", tz->prepro_memory_consump);
	printf ("Result of Thorup-Zwick dist(%d, %d) = %d\n", u, v, tz->dist);
	printf ("Time spent on dist Thorup-Zwick: %f sec\n", tz->dist_time);
	printf ("Query algorithm is executed %d times\n", tz->query_times);
	printf ("Memory usage of dist (bunch size) = %d KB\n", tz->dist_memory_consump);

	begin = clock();
	struct ssp_res *test = astar (graph, pp, u-offset, v-offset);
	end = clock();
	test->dist = test->S_f[v-offset].sp_est;
	printf ("Result of A*(%d,%d)=%d in time %f sec\n", u, v, test->dist, ((double)(end - begin) / CLOCKS_PER_SEC));
	return tz;
}

/**
 * run_tz - wrapper function for running Thorup-Zwick
 * @graph: graph with vertices and edges
 * @k: k integer
 * @u: source vertex u
 * @v: target vertex v
 * Calls Thorup-Zwick algorithm, and measures the spent RAM and CPU time
 * First it calls prepro, the preprocessing algorithm and then calls dist,
 * the query algorithm
 */
struct tz_res *run_tz (struct graph *graph, int k, int u, int v, int n, int m, int query_times)
{
	struct tz_res *tz = malloc (sizeof (struct tz_res));
	struct prepro *pp = malloc (sizeof (struct prepro));
	clock_t begin, end;

	tz->dist = 0, tz->dist_time = 0.0;
	tz->query_times = query_times;

	begin = clock();
	pp->success = false;
	while (!pp->success) {
	  pp = prepro (graph, k);
	}
	end = clock();
	tz->prepro_time = (double)(end - begin) / CLOCKS_PER_SEC;
	tz->prepro_memory_consump = get_vm_peak();

	for (int i = 0; i < tz->query_times; i++) {
		begin = clock();
		tz->dist += dist (&pp->nodes[u-offset], &pp->nodes[v-offset], pp->bunchlist);
		end = clock();
		tz->dist_time += (double)(end - begin) / CLOCKS_PER_SEC;
	}

	tz->dist = tz->dist / tz->query_times;
	tz->dist_time = tz->dist_time / tz->query_times;
	tz->dist_memory_consump += pp->bunchlist->bunch_size / 1000;
	tz->k = k;

	printf ("Time spent on prepro k=%d Thorup-Zwick: %f\n", tz->k, tz->prepro_time);
	printf ("vertices n=%d, edges m=%d\n", n, m);
	printf ("Memory usage of prepro = %d KB\n", tz->prepro_memory_consump);
	printf ("Result of Thorup-Zwick dist(%d, %d) = %d\n", u, v, tz->dist);
	printf ("Time spent on dist Thorup-Zwick: %f sec\n", tz->dist_time);
	printf ("Query algorithm is executed %d times\n", tz->query_times);
	printf ("Memory usage of dist (bunch size) = %d KB\n", tz->dist_memory_consump);

	begin = clock();
	struct ssp_res *test = astar (graph, pp, u-offset, v-offset);
	end = clock();
	test->dist = test->S_f[v-offset].sp_est;
	printf ("Result of A*(%d,%d)=%d in time %f sec\n", u, v, test->dist, ((double)(end - begin) / CLOCKS_PER_SEC));
	return tz;
}

/**
 * run_tz - wrapper function for running Thorup-Zwick
 * @graph: graph with vertices and edges
 * @k: k integer
 * @u: source vertex u
 * @v: target vertex v
 * Calls Thorup-Zwick algorithm, and measures the spent RAM and CPU time
 * First it calls prepro, the preprocessing algorithm and then calls dist,
 * the query algorithm
 */
struct tz_res *run_tz (struct graph *graph, int k, int u, int v, int n, int m, int query_times)
{
	struct tz_res *tz = malloc (sizeof (struct tz_res));
	struct prepro *pp = malloc (sizeof (struct prepro));
	clock_t begin, end;

	tz->dist = 0, tz->dist_time = 0.0;
	tz->query_times = query_times;

	begin = clock();
	pp->success = false;
	while (!pp->success) {
	  pp = prepro (graph, k);
	}
	end = clock();
	tz->prepro_time = (double)(end - begin) / CLOCKS_PER_SEC;
	tz->prepro_memory_consump = get_vm_peak();

	for (int i = 0; i < tz->query_times; i++) {
		begin = clock();
		tz->dist += dist (&pp->nodes[u-offset], &pp->nodes[v-offset], pp->bunchlist);
		end = clock();
		tz->dist_time += (double)(end - begin) / CLOCKS_PER_SEC;
	}

	tz->dist = tz->dist / tz->query_times;
	tz->dist_time = tz->dist_time / tz->query_times;
	tz->dist_memory_consump += pp->bunchlist->bunch_size / 1000;
	tz->k = k;

	printf ("Time spent on prepro k=%d Thorup-Zwick: %f\n", tz->k, tz->prepro_time);
	printf ("vertices n=%d, edges m=%d\n", n, m);
	printf ("Memory usage of prepro = %d KB\n", tz->prepro_memory_consump);
	printf ("Result of Thorup-Zwick dist(%d, %d) = %d\n", u, v, tz->dist);
	printf ("Time spent on dist Thorup-Zwick: %f sec\n", tz->dist_time);
	printf ("Query algorithm is executed %d times\n", tz->query_times);
	printf ("Memory usage of dist (bunch size) = %d KB\n", tz->dist_memory_consump);

	begin = clock();
	struct ssp_res *test = astar (graph, pp, u-offset, v-offset);
	end = clock();
	test->dist = test->S_f[v-offset].sp_est;
	printf ("Result of A*(%d,%d)=%d in time %f sec\n", u, v, test->dist, ((double)(end - begin) / CLOCKS_PER_SEC));
	return tz;
}

int main (int argc, char *argv[])
{

	if (argc == 2 && strcmp ("--help", argv[1]) == 0) {
		help ();
	} else if ((argc-1) < MIN_REQUIRED || (((argc-1) % MIN_REQUIRED) != 0)) {
		printf ("No input and output arguments or input/output does not match!\n");
		printf ("Number of arguments is %d\n\n", argc);
		help ();
		return EXIT_FAILURE;
	} else {
		for (int i = 1; i < argc; i += MIN_REQUIRED) {
			const char *fname_read = argv[i+1];
			offset = read_offset_in_file (fname_read);
			const char *fname_write = argv[i+2];
			const int u = atoi(argv[i+4]);
			const int v = atoi(argv[i+5]);
			const int query_times = atoi(argv[i+6]);
			struct graph_data *gd = read_vertices_and_edges (fname_read);
			if (u > gd->n || v > gd->n) {
				printf ("Source vertex u or/and target vertex v is/are invalid\n");
				printf ("Please consider a valid vertex that is <= n\n\n");
				help ();
				return EXIT_FAILURE;
			}
			struct graph *graph = init_graph (gd->n);
			read_from_file (graph, fname_read, offset);
			if (strcmp ("tz", argv[i]) == 0) {
				const int k = atoi(argv[i+3]);
				if (k <= 0) {
					printf ("k is <= 0! Please set k >= 1\n");
					help ();
					return EXIT_FAILURE;
				}
				struct tz_res *tz = run_tz (graph, k, u, v, gd->n, gd->m, query_times);
				write_to_csv (fname_write, fname_read, gd->n, gd->m, u, v,
							  tz, NULL, NULL, NULL);
			} else if (strcmp ("dj", argv[i]) == 0) {
				struct ssp_res *dijkstra = run_dijkstra (graph, u, v, gd->n, gd->m);
				write_to_csv (fname_write, fname_read, gd->n, gd->m, u, v,
							  NULL, dijkstra, NULL, NULL);
			} else if (strcmp ("djopt", argv[i]) == 0) {
				struct ssp_res *djopt = run_opt_dijkstra (graph, u, v, gd->n, gd->m);
				write_to_csv (fname_write, fname_read, gd->n, gd->m, u, v,
							  NULL, NULL, djopt, NULL);
			} else if (strcmp ("bdj", argv[i]) == 0) {
				struct ssp_res *bdj = run_bdj (graph, u, v, gd->n, gd->m);
				write_to_csv (fname_write, fname_read, gd->n, gd->m, u, v,
							  NULL, NULL, NULL, bdj);
			} else {
				printf ("The algorithm is not known!\n\n");
				help ();
			}
		};
	}
	return EXIT_SUCCESS;
}
