#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <map>
using namespace std;

#define DEBUG 0

struct edge{
	unsigned s;
	unsigned t;
	edge() = default;
	edge(unsigned a, unsigned b):s(a), t(b) {}
} ;
vector<vector<int>> adj;
typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges
	// vector<edge> edges;
	unsigned *ns;//ns[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *rd;
	unsigned *cd;//cumulative degree: (starts with 0) length=n+1
	unsigned *rcd;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj;//truncated list of neighbors
	unsigned *radj;//truncated list of neighbors
	unsigned *rank;//ranking of the nodes according to degeneracy ordering
	//unsigned *map;//oldID newID correspondance

	unsigned char *lab;//lab[i] label of node i
	unsigned **sub;//sub[l]: nodes in G_l

} specialsparse;

typedef struct {
	unsigned id;
	unsigned degree;
} iddegree;

iddegree *ig;
specialsparse *subg;
int K, *color, *ind, *loc, *C;
unsigned *cd0, *adj0, *dsub, *Index;

int cmp(const void* a, const void* b)
{
	iddegree *x = (iddegree*)a, *y = (iddegree*)b;
	return y->degree - x->degree;
}

int cmpadj(const void* a, const void* b)
{
	int *x = (int*)a, *y = (int*)b;
	return color[Index[*y]] - color[Index[*x]];
}

void freespecialsparse(specialsparse *g, unsigned char k) {
	unsigned char i;
	free(g->ns);
	for (i = 2; i < k + 1; i++) {
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
	free(g->edges);
	free(g->lab);
	free(g->cd);
	free(g->adj);
	free(g);
}

//Compute the maximum of three unsigned integers.
unsigned int max3(unsigned int a, unsigned int b, unsigned int c);
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}
int* v2m;
bool* st;
int* visCNT;
int* Candidate;
int CandidateCnt;
int partiteSize;
int* partiteVertexCSR;
int* partiteVertexSplit;
vector<int> partiteSetSize;

int N = 4, KClique;

vector<int> res;
size_t dfs_count;
int cliNum;
int* vMap;
int vMapCnt;
vector<vector<int>> allRes;

int* cnt;
// bool checkSt[4];
// int checkVertex[22470];
map<int, vector<vector<int>>> his;
map<int, bool> st2;
int dfs_cnt;
#include <chrono>
#include <fstream>
#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast< \
std::chrono::microseconds>(Get_Time() - start).count() / (float)1000
#define Print_Time(str, start) std::cout << str << Duration(start) << " ms" << std::endl
map<string, size_t> get_index_mem() {
	FILE* fp = fopen("/proc/self/status", "r");
	char line[128];
	map<string, size_t> res;
	while (fgets(line, 128, fp) != NULL)
	{
		//        if (strncmp(line, "VmPeak", 2) == 0)
		//        {
		//            cout << line << endl;
		////            printf("当前进程占用虚拟内存大小为：%d KB\n", atoi(line + 6));
		//        }
		if (strncmp(line, "VmRSS:", 6) == 0) {
			string p = line;
			res["now"] = size_t(stoull(p.substr(6)));
			//            cout << line;
		}
		if (strncmp(line, "VmPeak:", 7) == 0) {
			string p = line;
			res["pk"] = size_t(stoull(p.substr(7)));
			//            cout << line;
		}
	}
	fclose(fp);
	return res;
}

specialsparse* readedgelist(char* edgelist) {
	specialsparse *g = (specialsparse *)malloc(sizeof(specialsparse));
	FILE *file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist, "r");
	fscanf(file, "%u %u", &(g->n), &(g->e));
	cout << g->n << " " << g->e << endl;
	g->edges = (edge *)malloc(g->e * sizeof(edge));

	char c;
	v2m = new int[g->n];
	st = new bool[g->n];
	visCNT = new int[g->n];
	Candidate = new int[g->n];
	partiteVertexCSR = new int[g->n];

	memset(visCNT, 0, sizeof(int) * g->n);
	memset(st, 0, sizeof(bool) * g->n);
	memset(Candidate, 0, sizeof(int) * g->n);

	map<int, int> v2partiteSet;
	partiteSetSize.resize(partiteSize);
	for (int i = 0; i < g->n; ++i) {
		int x, y;
		// fscanf(file, " %c %u %u", &c, &x, &y);
		fscanf(file, "%u", &x);
		v2m[i] = x;
		// v2partiteSet[x].emplace_back(i);
		partiteSetSize[x]++;
	}
	// partiteVertexSplit = new int[v2partiteSet.size() + 1];
	// partiteVertexSplit[0] = 0;
	// for (int _ = 0; _ < v2partiteSet.size(); ++_)
	// 	partiteVertexSplit[_ + 1] = partiteVertexSplit[_] + v2partiteSet[_].size();
	// for (int _ = 0; _ < v2partiteSet.size(); ++_)
	// {
	// 	for (int j = partiteVertexSplit[_], _2 = 0; j < partiteVertexSplit[_ + 1]; ++j)
	// 	{
	// 		partiteVertexCSR[j] = v2partiteSet[_][_2++];
	// 	}
	// 	sort(partiteVertexCSR + partiteVertexSplit[_], partiteVertexCSR + partiteVertexSplit[_ + 1]);
	// }

	// vector<int> R2(g->n);
	// vector<int> R3(g->n);
	// string R2Path = edgelist;
	// string R3Path = edgelist;
	// R2Path += "R2";
	// R3Path += "R3";
	// FILE * file2 = fopen(R2Path.c_str(), "r");
	// for (int i = 0; i < g->n; ++i)
	// {
	// 	fscanf(file2, "%u", &(R2[i]));
	// }
	// fclose(file2);
	//
	// file2 = fopen(R3Path.c_str(), "r");
	// for (int i = 0; i < g->n; ++i)
	// {
	// 	fscanf(file2, "%u", &(R3[i]));
	// }
	// fclose(file2);


	adj.resize(g->n);

	for (int i = 0; i < g->e; ++i)
	{
		int s, t;
		fscanf(file, "%u %u", &s, &t);
		adj[s].emplace_back(t);
		adj[t].emplace_back(s);
	}
	fclose(file);

	vector<int> goodV(g->n);
	for (int i = 0; i < g->n; ++i)
	{
		map<int, bool> goodSt;
		for (auto j : adj[i])
		{
			goodSt[v2m[j]] = true;
		}
		if (goodSt.size() != partiteSize - 1)
			adj[i].clear();
	}

	file = fopen(edgelist, "r");
	fscanf(file, "%u %u", &(g->n), &(g->e));
	for (int i = 0; i < g->n; ++i) {
		int x;
		fscanf(file, "%u", &x);
	}
	int cnt_ = 0;
	for (int i = 0; i < g->e; ++i)
	{
		int s, t;
		fscanf(file, "%u %u", &s, &t);
		// if (!R2[s] || !R2[t]) continue;
		// if (adj[s].empty() || adj[t].empty()) continue;
		g->edges[cnt_].s = s;
		g->edges[cnt_].t = t;
		cnt_++;
	}
//	g->n++;
	g->e = cnt_;
	cout << g->e << endl;
//	g->edges = realloc(g->edges, g->e * sizeof(edge));

	return g;
}

void relabel(specialsparse *g) {
	unsigned i, source, target, tmp;
	vector<int> new_v2m(g->n);
	for (i = 0; i < g->n; i++)
	{
		new_v2m[g->rank[i]]  = v2m[i];
	}
	for (i = 0; i < g->n; i++)
	{
		v2m[i] = new_v2m[i];
	}
	for (i = 0; i < g->e; i++) {
		source = g->rank[g->edges[i].s];
		target = g->rank[g->edges[i].t];
		if (source < target) {
			tmp = source;
			source = target;
			target = tmp;
			//cout << "SWAP " << g->edges[i].s << " " << source << " " << g->edges[i].t << " " << target << endl;
			//swap(v2m[source], v2m[target]);
		}
		//swap(v2m[g->edges[i].s], v2m[source]);
		//swap(v2m[g->edges[i].t], v2m[target]);
		g->edges[i].s = source;
		g->edges[i].t = target;
	}
}

///// CORE ordering /////////////////////

typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(unsigned n_max) {
	unsigned i;
	bheap *heap = (bheap* )malloc(sizeof(bheap));

	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = (unsigned* )malloc(n_max * sizeof(unsigned));
	for (i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = (keyvalue* )malloc(n_max * sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap, unsigned i, unsigned j) {
	keyvalue kv_tmp = heap->kv[i];
	unsigned pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

void bubble_up(bheap *heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->kv[j].value > heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	unsigned i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n) {
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap, keyvalue kv) {
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_up(heap, heap->n - 1);
}

void update(bheap *heap, unsigned key) {
	unsigned i = heap->pt[key];
	if (i != -1) {
		((heap->kv[i]).value)--;
		bubble_up(heap, i);
	}
}

keyvalue popmin(bheap *heap) {
	keyvalue min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_down(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(unsigned n, unsigned *v) {
	unsigned i;
	keyvalue kv;
	bheap* heap = construct(n);
	for (i = 0; i < n; i++) {
		kv.key = i;
		kv.value = v[i];
		insert(heap, kv);
	}
	return heap;
}

void freeheap(bheap *heap) {
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void ord_core(specialsparse* g) {
	unsigned i, j, r = 0, n = g->n;
	keyvalue kv;
	bheap *heap;

	unsigned *d0 = (unsigned* )calloc(g->n, sizeof(unsigned));
	unsigned *cd0 = (unsigned* )malloc((g->n + 1) * sizeof(unsigned));
	unsigned *adj0 = (unsigned* )malloc(2 * g->e * sizeof(unsigned));
	for (i = 0; i < g->e; i++) {
		d0[g->edges[i].s]++;
		d0[g->edges[i].t]++;
	}
	cd0[0] = 0;
	for (i = 1; i < g->n + 1; i++) {
		cd0[i] = cd0[i - 1] + d0[i - 1];
		d0[i - 1] = 0;
	}
	for (i = 0; i < g->e; i++) {
		adj0[cd0[g->edges[i].s] + d0[g->edges[i].s]++] = g->edges[i].t;
		adj0[cd0[g->edges[i].t] + d0[g->edges[i].t]++] = g->edges[i].s;
	}

	heap = mkheap(n, d0);

	g->rank = (unsigned* )malloc(g->n * sizeof(unsigned));
	for (i = 0; i < g->n; i++) {
		kv = popmin(heap);
		g->rank[kv.key] = n - (++r);
		for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++) {
			update(heap, adj0[j]);
		}
	}
	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

//////////////////////////
//Building the special graph structure
void mkspecial(specialsparse *g, unsigned char k) {
	unsigned i, ns, max;
	unsigned *d, *sub;
	unsigned char *lab;

	d = (unsigned *)calloc(g->n, sizeof(unsigned));
	g->rd = (unsigned *)calloc(g->n, sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		d[g->edges[i].s]++;
		g->rd[g->edges[i].t]++;
	}

	g->cd = (unsigned* )malloc((g->n + 1) * sizeof(unsigned));
	g->rcd = (unsigned *)malloc((g->n + 1) * sizeof(unsigned));
	ns = 0;
	g->cd[0] = 0;
	g->rcd[0] = 0;
	max = 0;
	int max2 = 0;
	sub = (unsigned* )malloc(g->n * sizeof(unsigned));
	vector<pair<int, int>> sub_rk(g->n);
	lab = (unsigned char* )malloc(g->n * sizeof(unsigned char));
	for (i = 1; i < g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		g->rcd[i] = g->rcd[i - 1] + g->rd[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		max2 = (max2 > g->rd[i - 1]) ? max2 : g->rd[i - 1];
		// sub_rk[ns++] = {g->rank[i - 1], i - 1};
		sub[ns++] = i - 1;
		d[i - 1] = 0;
		g->rd[i - 1] = 0;
		lab[i - 1] = k;
	}
	// sort(sub_rk.begin(), sub_rk.end(), [](pair<int, int>& a, pair<int, int>& b)
	// {
	// 	return a.first < b.first;
	// });
	// for (i = 0; i < sub_rk.size(); i++)
	// {
	// 	sub[i] = sub_rk[i].second;
	// 	cout << sub_rk[i].second << endl;
	// }
	// printf("max degree = %u\n", max);

	subg = (specialsparse* )malloc(sizeof(specialsparse));
	subg->edges = (edge* )malloc((max*(max-1)/2) * sizeof(edge));
	subg->cd = (unsigned* )malloc((max + 1) * sizeof(unsigned));
	subg->adj = (unsigned* )malloc((max*(max - 1) / 2) * sizeof(unsigned));
	subg->ns = (unsigned* )malloc((k + 1) * sizeof(unsigned));

	subg->d = (unsigned** )malloc((k + 1) * sizeof(unsigned*));
	subg->sub = (unsigned** )malloc((k + 1) * sizeof(unsigned*));
	subg->lab = (unsigned char* )malloc(max * sizeof(unsigned char));

	for (int i = 2; i <= k; i++) {
		subg->d[i] = (unsigned* )malloc(max * sizeof(unsigned));
		subg->sub[i] = (unsigned* )malloc(max * sizeof(unsigned));
	}

	C = (int* )malloc(max * sizeof(int));
	ig = (iddegree* )malloc(max * sizeof(iddegree));
	adj0 = (unsigned* )malloc(2 * (max*(max-1)/2) * sizeof(unsigned));
	ind = (int* )malloc(g->n * sizeof(int));
	memset(ind, -1, g->n * sizeof(int));
	loc = (int* )malloc(max * sizeof(int));
	dsub = (unsigned *)calloc(max, sizeof(unsigned));
	Index = (unsigned* )malloc(max * sizeof(unsigned));
	cout << max << endl;
	color = (int* )malloc(max * sizeof(int));
	cd0 = (unsigned* )malloc((max + 1) * sizeof(unsigned));

	g->adj = (unsigned* )malloc(g->e * sizeof(unsigned));
	g->radj = (unsigned *)malloc(g->e * sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
		g->radj[g->rcd[g->edges[i].t] + g->rd[g->edges[i].t]++] = g->edges[i].s;
	}

	g->ns = (unsigned* )malloc((k + 1) * sizeof(unsigned));
	g->ns[k] = ns;

	g->d = (unsigned** )malloc((k + 1) * sizeof(unsigned*));
	g->sub = (unsigned** )malloc((k + 1) * sizeof(unsigned*));
	for (i = 2; i < k; i++) {
		g->d[i] = (unsigned* )malloc(g->n * sizeof(unsigned));
		g->sub[i] = (unsigned* )malloc(max * sizeof(unsigned));
	}
	g->d[k] = d;
	g->sub[k] = sub;

	g->lab = lab;
	map<int, bool> mps;

#if DEBUG
	for (int i = 0; i < g->n; i++)
	{
		cout << "#" << i << " ";
		auto end = g->cd[i] + g->d[partiteSize][i];
		auto rend = g->rcd[i] + g->rd[i];
		int ne;
		for (ne = g->cd[i]; ne < end; ne++)
			cout << g->adj[ne] << " ";
		cout << "| ";
		for (ne = g->rcd[i]; ne < rend; ne++)
			cout << g->radj[ne] << " ";
		cout << endl;
	}
#endif
}

void mkspecial_sub(specialsparse *g, unsigned char k) {
	unsigned i, ns, max;
	unsigned *d, *sub;
	unsigned char *lab;
	for (int i = 0; i < g->n; i++)
		g->d[k][i] = 0;
	for (i = 0; i < g->e; i++) {
		g->d[k][g->edges[i].s]++;

	}

	ns = 0;
	g->cd[0] = 0;
	max = 0;

	for (i = 1; i < g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + g->d[k][i - 1];
		max = (max > g->d[k][i - 1]) ? max : g->d[k][i - 1];
		g->sub[k][ns++] = i - 1;
		g->d[k][i - 1] = 0;
		g->lab[i - 1] = k;
	}

	for (i = 0; i < g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + g->d[k][g->edges[i].s]++] = g->edges[i].t;
	}

	g->ns[k] = ns;
}

int MaxColorNum;

void dfs(specialsparse *g, int u, int p, unsigned long long *n)
{
	if (u == K)
	{
#if DEBUG
		cout << "cliNum #" << ++cliNum << ": ";
		for (auto vv : res)
		{
			cout << vv << "(" << v2m[vv] << ") ";
		}
		cout << endl;

		// allRes.emplace_back(res);
		int sum = 0;
		for (auto vv : res) sum += vv;
		auto bk = res;
		sort(bk.begin(), bk.end());
		if (his[sum].empty()) his[sum].emplace_back(bk);
		else
		{
			for (auto vv : his[sum])
			{
				bool pass = true;
				for (int j = 0; j < vv.size(); j++)
				{
					if (vv[j] != bk[j])
					{
						pass = false;
						break;
					}
				}
				if (pass)
				{
					cout << "cliNum #" << *n << " same!!!" << endl;
					exit(0);
				}
			}
			his[sum].emplace_back(bk);
		}
#endif
		(*n)++;
		return;
	}
	int i, j, end, rend;
	for (i = p; i < CandidateCnt && CandidateCnt - i >= K - N - p; ++i)
	{
		bool flag = false;
		end = g->cd[Candidate[i]] + g->d[partiteSize][Candidate[i]];
		rend = g->rcd[Candidate[i]] + g->rd[Candidate[i]];

		for (j = g->cd[Candidate[i]]; j < end; j++) visCNT[g->adj[j]]++;
		for (j = g->rcd[Candidate[i]]; j < rend; j++) visCNT[g->radj[j]]++;
		if (u - cnt[v2m[Candidate[i]]] == visCNT[Candidate[i]])
		{
			flag = true;
			res.emplace_back(Candidate[i]);
			cnt[v2m[Candidate[i]]]++;
#if DEBUG
			// cout << "in " << dfs_cnt << " " << i << " add " << Candidate[i] << endl;
			dfs_cnt++;
#endif
			dfs(g, u + 1, i + 1, n);
#if DEBUG
			dfs_cnt--;
			// cout << "back " << dfs_cnt << " " << i << " del " << Candidate[i] << endl;
#endif

			cnt[v2m[Candidate[i]]]--;
		}
		for (j = g->cd[Candidate[i]]; j < end; j++) visCNT[g->adj[j]]--;
		for (j = g->rcd[Candidate[i]]; j < rend; j++) visCNT[g->radj[j]]--;
		if (flag) res.pop_back();
	}
}

void kclique(unsigned l, specialsparse *origin, specialsparse *g, unsigned long long *n) {
	unsigned i, j, k, end,end_, rend,rend_, u, v, w, ne;
	if (l == 2) {
		for (auto &_ : res) vMap[v2m[_]] = _;

		for (auto &_ : res) {
			end = origin->cd[_] + origin->d[partiteSize][_];
			rend = origin->rcd[_] + origin->rd[_];
			for (ne = origin->cd[_]; ne < end; ne++) visCNT[origin->adj[ne]]++;
			for (ne = origin->rcd[_]; ne < rend; ne++) visCNT[origin->radj[ne]]++;
		}

		for (i = 0; i < g->ns[2]; i++) {//list all edges
			u = g->sub[2][i];
			res.emplace_back(loc[u]);
			vMap[v2m[loc[u]]] = loc[u];
			// for (auto &ne : adj[u]) visCNT[ne]++;

			end = origin->cd[loc[u]] + origin->d[partiteSize][loc[u]];
			rend = origin->rcd[loc[u]] + origin->rd[loc[u]];
			for (ne = origin->cd[loc[u]]; ne < end; ne++) visCNT[origin->adj[ne]]++;
			for (ne = origin->rcd[loc[u]]; ne < rend; ne++) visCNT[origin->radj[ne]]++;

			end = g->cd[u] + g->d[2][u];

			for (j = g->cd[u]; j < end; j++) {
				v = loc[g->adj[j]];
				vMap[v2m[v]] = v;
				res.emplace_back(v);
				// cout << "!!!!!!!!!!!!!!!cliNum #" << ++cliNum << ": "; for (auto _ : res) cout << _ << " "; cout << endl;
				end_ = origin->cd[v] + origin->d[partiteSize][v];
				rend_ = origin->rcd[v] + origin->rd[v];
				for (ne = origin->cd[v]; ne < end_; ne++) visCNT[origin->adj[ne]]++;
				for (ne = origin->rcd[v]; ne < rend_; ne++) visCNT[origin->radj[ne]]++;

				for (int _ = 0; _ < origin->n; _++)
				{
					if (visCNT[_] == partiteSize - 1 && _ > vMap[v2m[_]]) // v是团中某三个点的邻居
					{
						Candidate[CandidateCnt++] = _;
#if DEBUG
						st2[Candidate[CandidateCnt - 1]] = true;
						cout << Candidate[CandidateCnt - 1] << "(" << v2m[Candidate[CandidateCnt - 1]] << ") ";
#endif
					}
				}
#if DEBUG
				cout << endl;
#endif

				for (int I = 0; I < partiteSize; I++) cnt[I] = 1;
				dfs(origin, partiteSize, 0, n);
				// (*n)++;
				CandidateCnt = 0;
				res.pop_back();
				end_ = origin->cd[v] + origin->d[partiteSize][v];
				rend_ = origin->rcd[v] + origin->rd[v];
				for (ne = origin->cd[v]; ne < end_; ne++) visCNT[origin->adj[ne]]--;
				for (ne = origin->rcd[v]; ne < rend_; ne++) visCNT[origin->radj[ne]]--;
			}

			res.pop_back();
			end = origin->cd[loc[u]] + origin->d[partiteSize][loc[u]];
			rend = origin->rcd[loc[u]] + origin->rd[loc[u]];
			for (ne = origin->cd[loc[u]]; ne < end; ne++) visCNT[origin->adj[ne]]--;
			for (ne = origin->rcd[loc[u]]; ne < rend; ne++) visCNT[origin->radj[ne]]--;
		}
		for (auto &_ : res) {
			end = origin->cd[_] + origin->d[partiteSize][_];
			rend = origin->rcd[_] + origin->rd[_];
			for (ne = origin->cd[_]; ne < end; ne++) visCNT[origin->adj[ne]] = 0;
			for (ne = origin->rcd[_]; ne < rend; ne++) visCNT[origin->radj[ne]] = 0;
		}
		return;
	}

	if (l == partiteSize)
	{
		for (i = 0; i < g->ns[l]; i++) {
			u = g->sub[l][i];

			g->ns[l - 1] = 0;
			end = g->cd[u] + g->d[l][u];
			for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
				v = g->adj[j];
				if (g->lab[v] == l) {		//equal to if(1)
					g->lab[v] = l - 1;
					g->sub[l - 1][g->ns[l - 1]++] = v;
					g->d[l - 1][v] = 0;//new degrees
				}
			}
			if (g->ns[l - 1] < 2)
			{
				for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
					v = g->sub[l - 1][j];
					g->lab[v] = l;
				}
				continue;
			}
			int cnt = -1, edge_num = 0;
			for (j = 0; j < g->ns[l - 1]; j++)
			{//reodering adjacency list and computing new degrees
				v = g->sub[l - 1][j];
				if (ind[v] == -1)
				{
					ind[v] = ++cnt;
					loc[cnt] = v;
					dsub[cnt] = 0;
				}
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					if (g->lab[w] == l - 1)
					{
						if (ind[w] == -1)
						{
							ind[w] = ++cnt;
							loc[cnt] = w;
							dsub[cnt] = 0;
						}
						edge_num++;
						dsub[ind[v]]++;
						dsub[ind[w]]++;
					}
				}
			}
			cd0[0] = 0;
			for (int i = 1; i < g->ns[l - 1] + 1; i++) {
				cd0[i] = cd0[i - 1] + dsub[i - 1];
				ig[i - 1].id = i - 1;
				ig[i - 1].degree = dsub[i - 1];
				dsub[i - 1] = 0;
			}

			for (j = 0; j < g->ns[l - 1]; j++)
			{
				color[j] = -1;
				v = g->sub[l - 1][j];
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					if (g->lab[w] == l - 1)
					{
						adj0[cd0[ind[v]] + dsub[ind[v]]++] = ind[w];
						adj0[cd0[ind[w]] + dsub[ind[w]]++] = ind[v];
					}
				}
			}

			qsort(ig, g->ns[l - 1], sizeof(ig[0]), cmp);

			for (int i = 0; i < g->ns[l - 1]; i++)
			{
				Index[ig[i].id] = i; // id -> i
				C[i] = 0;
			}

			color[0] = 0;
			int colorNum = 0;

			for (int i = 1; i < g->ns[l - 1]; i++)
			{
				int tmpdegree = ig[i].degree, tmpid = ig[i].id;

				for (int j = 0; j < tmpdegree; j++)
				{
					int now = Index[adj0[cd0[tmpid] + j]];
					if (color[now] != -1)
						C[color[now]] = 1;
				}
				for (int j = 0; j < ig[0].degree + 1; j++)
					if (C[j] == 0)
					{
						color[i] = j;
						colorNum = j > colorNum ? j : colorNum;
						break;
					}
				for (int j = 0; j < tmpdegree; j++)
				{
					int now = Index[adj0[cd0[tmpid] + j]];
					if (color[now] != -1)
						C[color[now]] = 0;
				}
			}

			int e_num = 0;
			for (j = 0; j < g->ns[l - 1]; j++)
			{
				v = g->sub[l - 1][j];
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					if (g->lab[w] == l - 1)
					{
						if (color[Index[ind[v]]] < color[Index[ind[w]]])
						{
							subg->edges[e_num].s = ind[w];
							subg->edges[e_num++].t = ind[v];
						}
						else if (color[Index[ind[v]]] == color[Index[ind[w]]])
						{
							if (ig[Index[ind[v]]].id < ig[Index[ind[w]]].id)
							{
								subg->edges[e_num].s = ind[v];
								subg->edges[e_num++].t = ind[w];
							}
							else
							{
								subg->edges[e_num].s = ind[w];
								subg->edges[e_num++].t = ind[v];
							}
						}
						else if (color[Index[ind[v]]] > color[Index[ind[w]]])
						{
							subg->edges[e_num].s = ind[v];
							subg->edges[e_num++].t = ind[w];
						}
					}
				}
			}

			subg->n = g->ns[l - 1];
			subg->e = edge_num;
			mkspecial_sub(subg, l - 1);

            if (colorNum + 1 > MaxColorNum) MaxColorNum = colorNum + 1;
			// printf("Color Nums: %d\n", colorNum + 1);
			res.emplace_back(u);

			kclique(l - 1, g, subg, n);
			res.pop_back();
			for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
				ind[loc[j]] = -1;
				v = g->sub[l - 1][j];
				g->lab[v] = l;
			}
		}
	}
	else
	{
		if (l > g->ns[l])
			return;
		for (int i = 0; i < g->ns[l]; i++) {
			u = g->sub[l][i];
			if (color[Index[u]] < l - 1)
				continue;
			g->ns[l - 1] = 0;
			end = g->cd[u] + g->d[l][u];
			for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
				v = g->adj[j];
				if (g->lab[v] == l) {
					g->lab[v] = l - 1;
					g->sub[l - 1][g->ns[l - 1]++] = v;
					g->d[l - 1][v] = 0;//new degrees
				}
			}
			for (j = 0; j < g->ns[l - 1]; j++) {//reodering adjacency list and computing new degrees
				v = g->sub[l - 1][j];
				end = g->cd[v] + g->d[l][v];
				int Index = g->cd[v];
				for (k = g->cd[v]; k < end; k++) {
					w = g->adj[k];
					if (g->lab[w] == l - 1) {
						g->d[l - 1][v]++;
					}
					else {
						g->adj[k--] = g->adj[--end];
						g->adj[end] = w;
					}

				}
			}
			res.emplace_back(loc[u]);
			kclique(l - 1, origin, g, n);
			res.pop_back();
			for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
				v = g->sub[l - 1][j];
				g->lab[v] = l;
			}
		}
	}
}
string getTime()
{
	time_t timep;
	time (&timep);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep));
	return tmp;
}
int main(int argc, char** argv) {
	cout << getTime() << endl;
	auto mp = get_index_mem();
	auto m1 = mp["pk"];
	specialsparse* g;
	partiteSize = atoi(argv[1]);
	unsigned char k = atoi(argv[2]);
	// unsigned char k = 7;
	K = k;
	unsigned long long n;
	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	printf("Reading edgelist from file %s\n", argv[3]);
	// g = readedgelist(argv[3]);
	// facebook amazon
	// char pa[100] = "../Dataset/amazon.csv";
	// char pa[100] = "../data/test";
	g = readedgelist(argv[3]);
	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);
	printf("Number of partiteSize = %u\n", partiteSize);
	cnt = new int[partiteSize];
	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Building the graph structure\n");

	ord_core(g);
	relabel(g);
	mkspecial(g, partiteSize);

	// for (int i = 0; i < g->n; ++i) {
	// 	cout << i << " " << Index[i] << " " << color[Index[i]] << endl;
	// }

	mp = get_index_mem();
	auto graphSize = mp["pk"];
	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	// t2 = time(NULL);
	//
	// printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	// t1 = t2;
	auto ti= Get_Time();
	printf("Iterate over all cliques\n");

	vMap = new int[g->n];
	memset(vMap, -1, g->n * sizeof(int));

	n = 0;
	kclique(partiteSize, g, g, &n);

	printf("Number of %u-cliques: %llu\n", k, n);
	Print_Time("All Time: ", ti);
	cout << "dfs_count: " << dfs_count << endl;
	mp = get_index_mem();
	auto m2 = mp["pk"];
	fstream fs = fstream("../base_res.csv", ios::app);
	fs << "comb," << (int)k << "," << n << "," << Duration(ti) << "," << m2 - m1 << "," << graphSize - m1 << endl;
	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

#if DEBUG
	for (auto& i : allRes) sort(i.begin(), i.end());
	sort(allRes.begin(), allRes.end());
	cout << endl;
	for (int i = 0; i < allRes.size(); i++)
	{
		cout << "#" << i + 1 << ": ";
		for (auto j : allRes[i])
		{
			cout << j << " ";
		}
		cout << endl;
	}
#endif
	freespecialsparse(g, partiteSize);

	free(color);
	free(Index);

	// printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}
