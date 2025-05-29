#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <ctime>
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
	unsigned value;
	unsigned degree;
} idrank;

int *color;
unsigned *Index;

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
string getTime()
{
	time_t timep;
	time (&timep);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep));
	return tmp;
}

int cmp_core_degree(const void* a, const void* b)
{
	idrank *x = (idrank*)a, *y = (idrank*)b;
	if (x->value != y->value)
		return y->value - x->value;
	else
		return y->degree - x->degree;
}

int cmp_color(const void* a, const void* b) {
  int u = *(int*)a, v = *(int*)b;
  if (color[Index[u]] < color[Index[v]]) return false;
  if (color[Index[u]] == color[Index[v]] && Index[u] > Index[v]) return false; // ir[Index[v]].id
  return true;
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
int* vec;
int vecSize;
int minSet;
map<int, int> partiteCount;
#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

specialsparse* readedgelist(char* edgelist) {
	specialsparse *g =(specialsparse*)malloc(sizeof(specialsparse));
	FILE *file;

	g->n = 0;
	g->e = 0;
	v2m = new int[g->n];

	file = fopen(edgelist, "r");
	fscanf(file, "%u", &g->n);
	fscanf(file, "%u", &g->e);
	cout << "Origin |V|:" << g->n << endl;
	cout << "Origin |E|:" << g->e << endl;
	vector<vector<int>> v2partiteSet(partiteSize);
	v2m = new int[g->n];
	for (int i = 0; i < g->n; i++) {
		fscanf(file, "%u", &v2m[i]);
		partiteCount[v2m[i]]++;
		// cout << v2m[i] << " " << partiteSize << endl;
		// assert(v2m[i] < partiteSize);
		v2partiteSet[v2m[i]].emplace_back(i);
	}

	cout << 1 << endl;
	for (auto& i : v2partiteSet)
	{
		i.shrink_to_fit();
	}
	int cnt = 0;
	unsigned e1 = NLINKS;
	g->edges = (edge*)malloc(e1 * sizeof(edge));
	for (int i = 0; i < g->e; i++)
	{
		int s, t;
		fscanf(file, "%u %u", &s, &t);
		g->edges[cnt++] = edge(s, t);

		if (cnt == e1) {
			e1 += NLINKS;
			g->edges = (edge*)realloc(g->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	cout << "After Del |E|:" << g->e << endl;

	for (auto& i : v2partiteSet)
	{
		for (int v1 = 0; v1 < i.size(); v1++)
			for (int v2 = v1 + 1; v2 < i.size(); v2++)
			{
				g->edges[cnt++] = edge(i[v1], i[v2]);
				if (cnt == e1) {
					e1 += NLINKS;
					g->edges = (edge*)realloc(g->edges, e1 * sizeof(edge));
				}
			}
	}
	g->e = cnt;
	g->edges = (edge*)realloc(g->edges, g->e * sizeof(edge));

	cout << "After Add |E|:" << g->e << endl;
	for (auto i : partiteCount)
	{
		cout << "#" << i.first << " " << i.second << endl;
	}
	return g;
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
	bheap *heap = (bheap *)malloc(sizeof(bheap));

	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = (unsigned *)malloc(n_max * sizeof(unsigned));
	for (i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = (keyvalue *)malloc(n_max * sizeof(keyvalue));
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

void update(bheap *heap, unsigned key,keyvalue kv) {
	unsigned i = heap->pt[key];
	if (i != -1) {
		if(kv.value < (heap->kv[i]).value)
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
void ord_color_relabel(specialsparse* g) {
	unsigned i, j, r = 0, N = g->n,maxdegree = 0;
	keyvalue kv;
	bheap *heap;

	idrank *ir = (idrank *)malloc(g->n * sizeof(idrank));
	unsigned *d0 = (unsigned *)calloc(g->n, sizeof(unsigned));
	unsigned *cd0 = (unsigned *)malloc((g->n + 1) * sizeof(unsigned));
	unsigned *adj0 = (unsigned *)malloc(2 * g->e * sizeof(unsigned));
	for (i = 0; i < g->e; i++) {
		d0[g->edges[i].s]++;
		d0[g->edges[i].t]++;
	}
	cd0[0] = 0;
	for (i = 1; i < g->n + 1; i++) {
		cd0[i] = cd0[i - 1] + d0[i - 1];

		maxdegree = (d0[i - 1] > maxdegree) ? d0[i - 1] : maxdegree;
		d0[i - 1] = 0;
	}
	for (i = 0; i < g->e; i++) {
		adj0[cd0[g->edges[i].s] + d0[g->edges[i].s]++] = g->edges[i].t;
		adj0[cd0[g->edges[i].t] + d0[g->edges[i].t]++] = g->edges[i].s;
	}

	heap = mkheap(N, d0);

	Index = (unsigned *)malloc(N * sizeof(unsigned));
	g->rank = (unsigned *)malloc(g->n * sizeof(unsigned));
	for (i = 0; i < g->n; i++) {
		kv = popmin(heap);
		ir[N-i-1].id = kv.key;
		//ir[i].rank = N - (r + 1);
		ir[N - i - 1].value = kv.value;
		ir[N - i - 1].degree = d0[kv.key];
		//core[kv.key] = kv.value;
		//Index[ir[N - i - 1].id] = N - i - 1;
		g->rank[kv.key] = N - (++r);
		for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++) {
			update(heap, adj0[j], kv);
		}
	}


	qsort(ir,N,sizeof(ir[0]),cmp_core_degree);
	for (int i = 0; i < N; i++)
	{
		// printf("id = %d value = %d degree = %d\n", ir[i].id,ir[i].value,ir[i].degree);
		Index[ir[i].id] = i;
	}

	// printf("after -----------\n");
	//
	// for (int i = 0; i < N; i++)
	// {
	// 	printf("id = %d value = %d degree = %d\n", ir[i].id, ir[i].value, ir[i].degree);
	// 	//Index[ir[i].id] = i;
	// }


	//color ordering
	color = (int *)malloc(N * sizeof(int));
	memset(color, -1, sizeof(int)*N);

	int *C = (int*)malloc((maxdegree+1) * sizeof(int));
	memset(C, 0, sizeof(int)*(maxdegree + 1));
	color[0] = 0;
	int colorNum = 1;


	for (int i = 1; i < N; i++)
	{
		int tmpdegree = ir[i].degree, tmpid = ir[i].id;
		for (int j = 0; j < tmpdegree; j++)
		{
			int now = Index[adj0[cd0[tmpid] + j]];
			if (color[now] != -1)
				C[color[now]] = 1;
		}
		for (int j = 0; j < maxdegree + 1; j++)
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
	printf("color number = %d max degree = %d\n", colorNum,maxdegree);

	//relabel
	for (int i = 0; i < g->e; i++)
	{
		if (color[Index[g->edges[i].s]] < color[Index[g->edges[i].t]])
		{
			int tmp = g->edges[i].s;
			g->edges[i].s = g->edges[i].t;
			g->edges[i].t = tmp;
		}
		else if (color[Index[g->edges[i].s]] == color[Index[g->edges[i].t]])
		{
			if (ir[Index[g->edges[i].s]].id > ir[Index[g->edges[i].t]].id)
			{
				int tmp = g->edges[i].s;
				g->edges[i].s = g->edges[i].t;
				g->edges[i].t = tmp;
			}
		}

	}

  if(false) { // added to output the color ordering and analyze its dpm-value
    printf("Output color ordering in ord.DegenCol.txt");
    for (int i = 0; i < N; i++) d0[i] = i;
    qsort(d0,N,sizeof(d0[0]),cmp_color); // color ordering
    for (int i = 0; i < N; i++) cd0[d0[i]] = i; // color rank

    FILE *file2 = fopen("ord.DegenCol.txt", "w");
    for (int i = 0; i < N; i++) {
      if(cd0[i] >= N) {
        printf("Prrrrroblem %d: %u\n", i, cd0[i]);
        continue;
      }
      fprintf(file2, "%u\n", cd0[i]);
    }
    fclose(file2);
  }

	free(C);
	free(ir);
	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

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
	g->cd = (unsigned *)malloc((g->n + 1) * sizeof(unsigned));
	g->rcd = (unsigned *)malloc((g->n + 1) * sizeof(unsigned));
	ns = 0;
	g->cd[0] = 0;
	g->rcd[0] = 0;
	max = 0;
	int max2 = 0;
	sub = (unsigned *)malloc(g->n * sizeof(unsigned));
	lab = (unsigned char *)malloc(g->n * sizeof(unsigned char));
	for (i = 1; i < g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		g->rcd[i] = g->rcd[i - 1] + g->rd[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		max2 = (max2 > g->rd[i - 1]) ? max2 : g->rd[i - 1];
		sub[ns++] = i - 1;
		d[i - 1] = 0;
		g->rd[i - 1] = 0;
		lab[i - 1] = k;
	}
	printf("max degree = %u\n", max);

	g->adj = (unsigned *)malloc(g->e * sizeof(unsigned));
	g->radj = (unsigned *)malloc(g->e * sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
		g->radj[g->rcd[g->edges[i].t] + g->rd[g->edges[i].t]++] = g->edges[i].s;
	}
	free(g->edges);


	g->ns = (unsigned *)malloc((k + 1) * sizeof(unsigned));
	g->ns[k] = ns;

	printf("max: %d\n", max);
	g->d = (unsigned **)malloc((k + 1) * sizeof(unsigned*));
	g->sub = (unsigned **)malloc((k + 1) * sizeof(unsigned*));
	for (i = 0; i <= k; i++) {
		g->d[i] = (unsigned *)malloc(g->n * sizeof(unsigned));
		g->sub[i] = (unsigned *)malloc(max * sizeof(unsigned));
	}
	g->d[k] = d;
	g->sub[k] = sub;

	g->lab = lab;
}
using namespace std;

int N = 4, K;

vector<int> res;
size_t dfs_count;
int cliNum;
int* vMap;
int vMapCnt;
vector<vector<int>> allRes;

map<int, vector<vector<int>>> his;
map<int, bool> st2;
int dfs_cnt;
int A, B;

int* CNT;
int MaxCnt;
void kclique(unsigned l, specialsparse *g, unsigned long long *n) {
	unsigned i, j, k, end,end_, rend,rend_, u, v, w, ne;
	// if (l == 2) {
	// 	for (i = 0; i < g->ns[2]; i++) {//list all edges
	// 		u = g->sub[2][i];
	// 		if (CNT[v2m[u]] + 1 > MaxCnt) continue;
	// 		end = g->cd[u] + g->d[2][u];
	// 		CNT[v2m[u]]++;
	// 		// res.emplace_back(u);
	// 		for (j = g->cd[u]; j < end; j++) {
	// 			if (CNT[v2m[g->adj[j]]] + 1 > MaxCnt) continue;
	// 			CNT[v2m[g->adj[j]]]++;
				bool ok = true;

				for (int kd = 0; kd < partiteSize; kd++)
				{
					// if (CNT[kd] == 0)
					// {
					// 	ok = false;
					// 	break;
					// }
					// 2
					if (CNT[kd] < A)
					{
						ok = false;
						break;
					}
				}

				if (ok)
				{
					// int mx = -1, mn = 0x3f3f3f;
					// for (int kd = 0; kd < partiteSize; kd++)
					// {
					// 	mx = max(mx, CNT[kd]);
					// 	mn = min(mn, CNT[kd]);
					// 	if (mx - mn > A) {
					// 		ok = false;
					// 		break;
					// 	}
					// }
					// if (ok)
						(*n)++;//listing here!!!
				}

	// 			CNT[v2m[g->adj[j]]]--;
	// 		}
	// 		CNT[v2m[u]]--;
	// 		// res.pop_back();
	// 	}
	// 	return;
	// }

	// if (l > g->ns[l])
	// 	return;

	for (i = 0; i < g->ns[l]; i++) {
		u = g->sub[l][i];
		// 2
		if (CNT[v2m[u]] >= B) {
			continue;
		}
		// if (color[Index[u]] < l - 1)
		// 	continue;
		g->ns[l - 1] = 0;
		end = g->cd[u] + g->d[l][u];
		for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
			v = g->adj[j];
			if (g->lab[v] == l) {
				g->lab[v] = l - 1;
				g->sub[l - 1][g->ns[l - 1]++] = v;
				g->d[l - 1][v] = 0;//new degrees
				// cout << "l-1(" << l - 1 << ") " << v << endl;
			}
		}

		for (j = 0; j < g->ns[l - 1]; j++) {//reodering adjacency list and computing new degrees
			v = g->sub[l - 1][j];
			end = g->cd[v] + g->d[l][v];
			int Index = g->cd[v];// , tol = g->cd[v];
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
		// if (CNT[v2m[u]] + 1 <=MaxCnt)
		// {
			CNT[v2m[u]]++;
			// res.emplace_back(u);
			kclique(l - 1, g, n); // Gv = g[l-1]
			// res.pop_back();
			CNT[v2m[u]]--;
		// }
		for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
			v = g->sub[l - 1][j];
			g->lab[v] = l;
		}
	}
}

int main(int argc, char** argv) {
	cout << getTime() << endl;
	auto mp = get_index_mem();
	auto m1 = mp["pk"];
	specialsparse* g;
	partiteSize = atoi(argv[1]);
	// partiteSize = 4;
	unsigned char k = atoi(argv[2]);
	// unsigned char k = 8;
	K = k;
	unsigned long long n;
	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;
	MaxCnt = K - partiteSize + 1;
	printf("Reading edgelist from file %s\n", argv[3]);
	g = readedgelist(argv[3]);
	A = atoi(argv[4]);
	B = atoi(argv[5]);
	// A = 1;
	// B = 2;
	// facebook amazon
	// char pa[100] = "../../Dataset/KPC_Dataset/congressZJ4.csv";
	// char pa[100] = "../data/test";
	// g = readedgelist(pa);
	CNT = new int[partiteSize];
	memset(CNT, 0, partiteSize * sizeof(int));
	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);
	printf("Number of partiteSize = %u\n", partiteSize);
	printf("A=%d B=%d\n", A, B);
	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Building the graph structure\n");

	ord_color_relabel(g);

	mkspecial(g, k);

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
	kclique(k, g, &n);

	printf("Number of %u-cliques: %llu\n", k, n);
	Print_Time("All Time: ", ti);
	cout << "dfs_count: " << dfs_count << endl;
	mp = get_index_mem();
	auto m2 = mp["pk"];
	fstream fs = fstream("../base_res.csv", ios::app);
	fs << argv[3] << ",base2," << g->n << "," << g->e << "," << (int)k << "," << n << "," << Duration(ti) << "," << m2 - m1 << "," << graphSize - m1 << "," << A << "," << B << endl;
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
	freespecialsparse(g, k);

	free(color);
	free(Index);

	// printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}
