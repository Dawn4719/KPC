#include <assert.h>
#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <map>
#include <iostream>
#define bzero(a, b)             memset(a, 0, b)
using namespace std;

#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast<std::chrono::microseconds>(Get_Time() - start).count()/(float)1000
#define Print_Time(str, start) std::cout << str << Duration(start) << "ms" << std::endl

struct edge{
	unsigned s;
	unsigned t;
	edge() = default;
	edge(unsigned s_, unsigned t_):s(s_), t(t_) {}
} ;

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg;


typedef struct {

	unsigned n;//number of nodes
	unsigned e;//number of edges
	vector<edge> edges;//list of edges

	unsigned *ns;//ns[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *cd;//cumulative degree: (start with 0) length=n+1
	unsigned *adj;//list of neighbors with lower degree

	unsigned char *lab;//lab[i] label of node i
	unsigned **sub;//sub[l]: nodes in G_l

} specialsparse;


void freespecialsparse(specialsparse *g, unsigned char k) {
	unsigned char i;
	free(g->ns);
	for (i = 2; i<k + 1; i++) {
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
	free(g->cd);
	free(g->adj);
	free(g);
}

//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a>b) ? a : b;
	return (a>c) ? a : c;
}
vector<int> partiteSet;
map<int, int> partiteCount;
int K;
specialsparse* readedgelist(char* edgelist) {
	specialsparse *g =(specialsparse*)malloc(sizeof(specialsparse));
	FILE *file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist, "r");
	fscanf(file, "%u", &g->n);
	fscanf(file, "%u", &g->e);
	partiteSet.resize(g->n);

	for (int i = 0; i < g->n; i++) {
		fscanf(file, "%u", &partiteSet[i]);
		partiteCount[partiteSet[i]]++;
	}
	g->edges.reserve(g->e);
	fclose(file);

	file = fopen(edgelist, "r");
	fscanf(file, "%u", &g->n);
	fscanf(file, "%u", &g->e);
	vector<vector<int>> v2partiteSet(partiteCount.size());
	for (int i = 0; i < partiteCount.size(); ++i) { v2partiteSet[i].reserve(partiteCount[i]); }
	int tmp_;
	for (int i = 0; i < g->n; i++) {
		fscanf(file, "%u %u", &tmp_, &tmp_);
		v2partiteSet[partiteSet[i]].emplace_back(i);
	}
	int cnt = 0;
	for (int i = 0; i < g->e; i++)
	{
		int s, t;
		fscanf(file, "%u %u", &s, &t);
		g->edges.emplace_back(s, t);
		cnt++;
	}
	g->e = cnt;

	g->edges.shrink_to_fit();
	fclose(file);

	return g;
}

//Building the special graph structure
void mkspecial(specialsparse *g, unsigned char k) {
	unsigned i, ns, max;
	unsigned *d, *sub;
	unsigned char *lab;

	d = (unsigned*)calloc(g->n, sizeof(unsigned));

	for (i = 0; i<g->e; i++) {
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}

	g->cd = (unsigned*)malloc((g->n + 1) * sizeof(unsigned));
	ns = 0;
	g->cd[0] = 0;
	max = 0;
	sub = (unsigned*)malloc(g->n * sizeof(unsigned));
	lab = (unsigned char*)malloc(g->n * sizeof(unsigned char));
	for (i = 1; i<g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		if (d[i - 1]>0) {
			max = (max>d[i - 1]) ? max : d[i - 1];
			sub[ns] = i - 1;
			ns++;
			d[i - 1] = 0;
			lab[i - 1] = k;
		}
	}

	g->adj = (unsigned*)malloc(2 * g->e * sizeof(unsigned));

	for (i = 0; i<g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
		g->adj[g->cd[g->edges[i].t] + d[g->edges[i].t]++] = g->edges[i].s;
	}
	// free(g->edges);

	g->ns = (unsigned*)malloc((k + 1) * sizeof(unsigned));
	g->ns[k] = ns;

	g->d = (unsigned**)malloc((k + 1) * sizeof(unsigned*));
	g->sub = (unsigned**)malloc((k + 1) * sizeof(unsigned*));
	for (i = 2; i<k; i++) {
		g->d[i] = (unsigned*)malloc(g->n * sizeof(unsigned));
		g->sub[i] = (unsigned*)malloc(max * sizeof(unsigned));
	}
	g->d[k] = d;
	g->sub[k] = sub;

	g->lab = lab;
}

/*For futur use in qsort_r.
int compare_r(void const *p_i, void const *p_j,void *p_deg){
unsigned i = *((unsigned*)p_i);
unsigned j = *((unsigned*)p_j);
unsigned *deg = (unsigned*)p_deg;
return (deg[i] < deg[j]) ? 1 : -1;
}
*/

void arg_bucket_sort(unsigned *key, unsigned n, unsigned *val) {
	unsigned i, j;
	static unsigned *c = NULL, *cc = NULL, *key2 = NULL;
	if (c == NULL) {
		c = (unsigned*)malloc(n * sizeof(unsigned));//count
		cc = (unsigned*)malloc(n * sizeof(unsigned));//cummulative count
		key2 = (unsigned*)malloc(n * sizeof(unsigned));//sorted array
	}
	bzero(c, n * sizeof(unsigned));

	for (i = 0; i<n; i++) {
		(c[val[key[i]]])++;
	}
	cc[0] = 0;
	for (i = 1; i<n; i++) {
		cc[i] = cc[i - 1] + c[i - 1];
		c[i - 1] = 0;
	}
	c[i - 1] = 0;

	for (i = 0; i<n; i++) {
		j = val[key[i]];
		key2[cc[j] + c[j]++] = key[i];
	}

}
vector<int> Clique;
vector<int> trianglesCnt;
map<int, bool> kinds;

void kclique(unsigned l, specialsparse *g, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;

	if (l == 2) {
		for (i = 0; i<g->ns[2]; i++) {//list all edges
			u = g->sub[2][i];
			end = g->cd[u] + g->d[2][u];
			for (j = g->cd[u]; j<end; j++) {
				v = g->adj[j];
				if (v<u) {
					for (auto _ : Clique)
					{
						kinds[partiteSet[_]] = true;
						// trianglesCnt[_]++;
					}
					kinds[partiteSet[u]] = true;
					kinds[partiteSet[v]] = true;
					// trianglesCnt[u]++;
					// trianglesCnt[v]++;

					cout << "Clique: ";
					for (auto _ : Clique)
					{
						cout << _ << " ";
					}
					cout << u << " " << v << endl;
					(*n)++;//listing here!!!

				}
			}
		}
		return;
	}

	arg_bucket_sort(g->sub[l], g->ns[l], g->d[l]);
	//qsort_r(g->sub[l],g->ns[l],sizeof(unsigned),compare_r,g->d[l]);//qsort and bucket sort leads to similar running time in practice

	for (i = 0; i<g->ns[l]; i++) {
		u = g->sub[l][i];
		Clique.emplace_back(u);
		g->ns[l - 1] = 0;
		end = g->cd[u] + g->d[l][u];
		for (j = g->cd[u]; j<end; j++) {//relabeling nodes and forming U'.
			v = g->adj[j];
			if (g->lab[v] == l) {
				g->lab[v] = l - 1;
				g->sub[l - 1][g->ns[l - 1]++] = v;
				g->d[l - 1][v] = 0;//new degrees
			}
		}
		for (j = 0; j<g->ns[l - 1]; j++) {//reodering adjacency list and computing new degrees
			v = g->sub[l - 1][j];
			end = g->cd[v] + g->d[l][v];
			for (k = g->cd[v]; k<end; k++) {
				w = g->adj[k];
				if (g->lab[w] == l - 1) {
					g->d[l - 1][w]++;
				}
				else {
					g->adj[k--] = g->adj[--end];
					g->adj[end] = w;
				}
			}
		}

		kclique(l - 1, g, n);
		Clique.pop_back();
		for (j = 0; j<g->ns[l - 1]; j++) {//moving u to last position in each entry of the adjacency list
			v = g->sub[l - 1][j];
			g->lab[v] = l;
			end = g->cd[v] + g->d[l - 1][v];
			for (k = g->cd[v]; k<end; k++) {
				w = g->adj[k];
				if (w == u) {
					g->adj[k] = g->adj[--end];
					g->adj[end] = w;
					g->d[l - 1][v]--;
					break;
				}
			}
		}
		// g->lab[u] = l + 1;
		if (l == K)
		{
			for (auto& _ : kinds)
				assert(_.second == true);
			if (kinds.size() == partiteCount.size())
				trianglesCnt[u] = *n;
			else
				trianglesCnt[u] = 0;
			*n = 0;
			kinds.clear();
		}
	}
}

int main(int argc, char** argv) {
	specialsparse* g;

	unsigned char k = 3;
	K = k;
	cout << "# K = " << K << endl;
	unsigned long long n;

	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	char dataset[100] = "../../Dataset/facebook.csv";

	g = readedgelist(dataset);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);
	printf("Number of partiteCount = %u\n", partiteCount.size());

	printf("Building the graph structure\n");

	mkspecial(g, k);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Iterate over all cliques\n");

	//qsort_r(g->sub[k],g->ns[k],sizeof(unsigned),compare_r,g->d[k]);
	//arg_bucket_sort(g->sub[k], g->ns[k], g->d[k]);
	n = 0;
	trianglesCnt.resize(g->n);
	kclique(k, g, &n);

	string R2 = dataset;
	R2 += "R2";
	fstream f (R2, fstream::out);
	for (int i = 0; i < g->n; ++i)
		f << trianglesCnt[i] << endl;
	f.close();

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	freespecialsparse(g, k);

	printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}