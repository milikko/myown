#include "thorupzwick.h"

/**
 * dist - distance query of (u, v) with the data structures
 * @u: vertex u
 * @v: vertex v
 * @bunchlist: bunches of all nodes
 * Running time: O(k) = O(1)
 * u and v are the two vertices whose distance is to be estimated.
 * If w in B(v) it returns d(w, u) + d(w, v). The distance d(w, u) = d(p_i(u), u) is
 * read directly from the data structure constructed in prepro. Also, d(w, v) = d(v, w)
 * at most k accesses as w in A_{k-1} and A_{k-1} \subseteq for every v in V
 */
int dist (struct node *u, struct node *v, struct bunchlist *bunchlist)
{
	struct node	*w;
	int i = 0;
	int result = 0;

	// w = u = p_0(u)
	w = u;

	// while w not in B(v)
	while (1) {
		struct node *out, *tmp;
		// checking if w in B(v), where .nodes rep. all vertices w in V
		HASH_FIND_INT(bunchlist->bunches[v->v_id].nodes, &w->v_id, out);
		if (out) {
			result += out->sp_est;
			break;
		} else {
			// i <- i + 1
			i += 1;
			// (u, v) <- (v, u)
			tmp = u;
			u = v;
			v = tmp;
			// w <- p_i (u)
			w = &bunchlist->bunches[u->v_id].piv[i];
			result = w->sp_est;
		}
	}

	return result;
}
