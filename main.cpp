//
// Created by ssunah on 12/3/18.
//

#include <algorithm>
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

int* partiteSet;
map<int, int> partiteCount;
int K;
vector<int> degree;
int interpartiteEdges;
vector<vector<int>> adj;
int N;
specialsparse* readedgelist(char* edgelist) {
	specialsparse *g =(specialsparse*)malloc(sizeof(specialsparse));
	FILE *file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist, "r");
	fscanf(file, "%u", &g->n);
	fscanf(file, "%u", &g->e);
	partiteSet = new int[g->n];

	for (int i = 0; i < g->n; i++) {
		fscanf(file, "%u", &partiteSet[i]);
		partiteCount[partiteSet[i]] = 0;
	}
	N = partiteCount.size();
	adj.resize(g->n);
	for (int i = 0; i < g->e; i++)
	{
		int s, t;
		fscanf(file, "%u %u", &s, &t);
		adj[s].emplace_back(t);
		adj[t].emplace_back(s);
	}
	fclose(file);

	map<int, int> kinds;
	int cnt_ = 0;
	// for (auto & i : adj)
	// {
	// 	for (auto j : i)
	// 	{
	// 		kinds[partiteSet[j]] = true;
	// 	}
	// 	if (kinds.size() == N - 1)
	// 	{
	// 		kinds.clear();
	// 		continue;
	// 	}
	// 	i.clear();
	// 	kinds.clear();
	// 	cnt_++;
	// }

	vector<int> R2(g->n);
	vector<int> R3(g->n);
	string R2Str = edgelist;
	string R3Str = edgelist;
	fstream f(R2Str + "R2", ios::in);
	for (int i = 0; i < g->n; i++)
		f >> R2[i];
	f.close();

	f = fstream(R3Str + "R3", ios::in);
	for (int i = 0; i < g->n; i++)
	{
		f >> R3[i];
	}
	f.close();

	file = fopen(edgelist, "r");
	fscanf(file, "%u", &g->n);
	fscanf(file, "%u", &g->e);
	cout << "Origin |E|:" << g->e << endl;
	vector<vector<int>> v2partiteSet(partiteCount.size());
	for (int i = 0; i < N; ++i) { v2partiteSet[i].reserve(partiteCount[i]); }
	int tmp_;
	for (int i = 0; i < g->n; i++) {
		fscanf(file, "%u", &tmp_);
		v2partiteSet[partiteSet[i]].emplace_back(i);
	}
	int cnt = 0;
	for (int i = 0; i < g->e; i++)
	{
		int s, t;
		fscanf(file, "%u %u", &s, &t);
		if (adj[s].empty() || adj[t].empty()) continue;
		// if (R2[s] + 1 < N || R2[t] + 1 < partiteCount.size()) continue;
		// if (!R3[s] || !R3[t]) continue;
		g->edges.emplace_back(s, t);
		partiteCount[partiteSet[s]]++;
		partiteCount[partiteSet[t]]++;
		cnt++;
	}
	g->e = cnt;
	fclose(file);
	cout << "After Del |E|:" << g->e << endl;
	// interpartiteEdges = 0;
	// for (auto& _ : v2partiteSet) {
	//     for (int s = 0; s < _.size(); ++ s) {
	//         for (int t = s + 1; t < _.size(); ++ t) {
	//         	if (adj[s].empty() || adj[t].empty()) continue;
	//         	if (R2[s] + 1 < N || R2[t] + 1 < partiteCount.size()) continue;
	//         	if (!R3[s] || !R3[t]) continue;
	//         	g->edges.emplace_back(_[s], _[t]);
	//         	g->e++;
	//             interpartiteEdges ++;
	//         }
	//     }
	// }
	g->edges.shrink_to_fit();
	cout << "After Add |E|:" << g->e << endl;
	for (auto i : partiteCount)
	{
		cout << "#" << i.first << " " << i.second << endl;
	}
	return g;
}
int beginPartite;
//Building the special graph structure
/*
void mkspecial(specialsparse *g, unsigned char k) {
	unsigned i, ns, max;
	unsigned *d, *sub;
	unsigned char *lab;

	d = (unsigned*)calloc(g->n, sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}

	// vector<pair<int, int>> partiteDegree(partiteCount.size());
	// for (i = 0; i < partiteDegree.size(); i++)
	// 	partiteDegree[i].second = i;
	// for (i = 0; i < g->n; ++i)
	// {
	// 	partiteDegree[partiteSet[i]].first += d[i];
	// }
	//
	// sort(partiteDegree.begin(), partiteDegree.end());
	//
	// for (i = g->e - interpartiteEdges; i < g->e; ++i)
	// {
	// 	d[g->edges[i].s]++;
	// 	d[g->edges[i].t]++;
	// }

	g->cd = (unsigned*)malloc((g->n + 1) * sizeof(unsigned));
	ns = 0;
	g->cd[0] = 0;
	max = 0;
	// sub = (unsigned*)malloc(g->n * sizeof(unsigned));
	sub = (unsigned*)malloc(8 * sizeof(unsigned));
	lab = (unsigned char*)malloc(g->n * sizeof(unsigned char));
	cout << endl;
	for (i = 1; i<g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		if (d[i - 1]>0) {
			max = (max>d[i - 1]) ? max : d[i - 1];
			// if (partiteSet[i - 1] == partiteDegree[0].second)
			// {
				sub[ns] = i - 1;
				ns++;
			// }
			d[i - 1] = 0;
			lab[i - 1] = k;
		}
	}
	// cout << "Max: " << max << endl;
	// cout << endl << "Minimal partite Set: " << partiteDegree[0].second << endl;
	// cout << "Minimal partite Set size: " << partiteCount[partiteDegree[0].second] << " " << ns << endl;
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
*/
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

vector<int> Clique;
map<int, int> kinds;
vector<int> visCNT;
vector<bool> st;
vector<int> C;
int Ccnt;

void dfs(int u, int p, unsigned long long *n)
{
	if (p + N > K)
	{
		(*n)++;
		return;
	}
	for (int i = p; i < Ccnt; ++i)
	{
		dfs(u + 1, i + 1, n);
	}
}

int calcu_cnt;
void arg_bucket_sort(unsigned *key, unsigned n, unsigned *val) {
	unsigned i, j;
	static unsigned *c = NULL, *cc = NULL, *key2 = NULL;
	if (c == NULL) {
		c = static_cast<unsigned*>(malloc(n * sizeof(unsigned)));//count
		cc = static_cast<unsigned*>(malloc(n * sizeof(unsigned)));//cummulative count
		key2 = static_cast<unsigned*>(malloc(n * sizeof(unsigned)));//sorted array
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

void kclique(unsigned l, specialsparse *g, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;

	if (l == 2) {
		for (auto& _ : Clique)
			for (auto& nei : adj[_])
			{
				if (st[nei] == false)
					visCNT[nei]++;
			}
		for (i = 0; i<g->ns[2]; i++) {//list all edges
			u = g->sub[2][i];
			end = g->cd[u] + g->d[2][u];
			for (auto &ne : adj[u]) {
				if (st[ne]) continue;
					visCNT[ne]++;
			}
			for (j = g->cd[u]; j<end; j++) {
				v = g->adj[j];
				if (v<u) {
					// cout << "Clique: ";
					// for (auto _ : Clique)
					// {
					// 	cout << _ << " ";
					// }
					// cout << u << " " << v << endl;
					for (auto &ne : adj[v]) {
						if (st[ne]) continue;
						visCNT[ne]++;
					}
					for (int _ = 0; _ < g->n; _++)
					{
						if (visCNT[_] >= N - 1)
						{
							C[Ccnt++] = _;
						}
					}
					(*n)++;
					dfs(0, 0, n);
					//listing here!!!
					for (auto &ne : adj[v]) {
						if (st[ne]) continue;
						visCNT[ne]--;
					}
				}
			}
			for (auto &ne : adj[u]) {
				if (st[ne]) continue;
				visCNT[ne]--;
			}
		}
		for (auto& _ : Clique)
		{
			for (auto& nei : adj[_])
			{
				visCNT[nei] = 0;
			}
		}
		return;
	}
	arg_bucket_sort(g->sub[l], g->ns[l], g->d[l]);
	for (i = 0; i<g->ns[l]; i++) {
		u = g->sub[l][i];
		if (l == 4)
		{
			st[u] = true;
			if (partiteSet[u] != 0) continue;
		}
		Clique.emplace_back(u);
		g->ns[l - 1] = 0;
		end = g->cd[u] + g->d[l][u];
		for (j = g->cd[u]; j<end; j++) {//relabeling nodes and forming U'.
			v = g->adj[j];
			// cout << v << endl;
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
		g->lab[u] = l + 1;
	}
}

int main(int argc, char** argv) {
	specialsparse* g;

	unsigned char k = 5;
	K = k;

	cout << "# K = " << K << endl;
	unsigned long long n;

	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	char dataset[100] = "../../Dataset/facebook.csv";
	// char dataset[100] = "../data/test";
	g = readedgelist(dataset);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);
	printf("Number of partiteCount = %u\n", partiteCount.size());

	printf("Building the graph structure\n");

	mkspecial(g, N);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Iterate over all cliques\n");

	//qsort_r(g->sub[k],g->ns[k],sizeof(unsigned),compare_r,g->d[k]);
	//arg_bucket_sort(g->sub[k], g->ns[k], g->d[k]);
	n = 0;

	for (int i = 0; i < N; i++) kinds[i] = false;
	visCNT.resize(g->n);
	C.resize(g->n);
	st.resize(g->n);
	kclique(N, g, &n);

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	freespecialsparse(g, N);

	printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}