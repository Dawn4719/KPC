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
#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

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
	edge *edges;//list of edges

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

int* v2m;
map<int, int> partiteCount;
int K;
int interpartiteEdges;
int partiteSize;
int MaxCnt;
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
		v2partiteSet[v2m[i]].emplace_back(i);
	}

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
				// cnt++;
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
int beginPartite;
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
	for (i = 0; i <= k; i++) {
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

int* CNT;
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
int A, B;
void kclique(unsigned l, specialsparse *g, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;

	if (l == 2) {
		for (i = 0; i<g->ns[2]; i++) {//list all edges
			u = g->sub[2][i];
			end = g->cd[u] + g->d[2][u];
			// if (CNT[v2m[u]] >= MaxCnt) continue;
			// 2
			// if (CNT[v2m[u]] >= B) {
			// 	continue;
			// }
			CNT[v2m[u]]++;

			for (j = g->cd[u]; j<end; j++) {
				v = g->adj[j];
				if (v<u) {
					if (CNT[v2m[v]] >= MaxCnt) continue;
	// 				// 2
	// 				// if (CNT[v2m[u]] >= B) {
	// 				// 	continue;
	// 				// }
					// CNT[v2m[v]]++;
					// Clique.emplace_back(v);
					// bool ok = true;
					//
					// for (int kd = 0; kd < partiteSize; kd++)
					// {
					// 	if (CNT[kd] == 0)
					// 	{
					// 		ok = false;
					// 		break;
					// 	}
					// 	// 2
					// 	// if (CNT[kd] < A)
					// 	// {
					// 	// 	ok = false;
					// 	// 	break;
					// 	// }
					// }

					// if (ok)
					// {
					// 	int mx = -1, mn = 0x3f3f3f;
					// 	for (int kd = 0; kd < partiteSize; kd++)
					// 	{
					// 		mx = max(mx, CNT[kd]);
					// 		mn = min(mn, CNT[kd]);
					// 		if (mx - mn > A) {
					// 			ok = false;
					// 			break;
					// 		}
					// 	}
					// 	if (ok)
							(*n)++;//listing here!!!
					// }

					// CNT[v2m[v]]--;
					// Clique.pop_back();
				}
			}
			CNT[v2m[u]]--;
		}
		return;
	}

	arg_bucket_sort(g->sub[l], g->ns[l], g->d[l]);
	//qsort_r(g->sub[l],g->ns[l],sizeof(unsigned),compare_r,g->d[l]);//qsort and bucket sort leads to similar running time in practice

	for (i = 0; i<g->ns[l]; i++) {
		u = g->sub[l][i];
		// 2
		// if (CNT[v2m[u]] >= B) {
		// 	continue;
		// }
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

		// if (CNT[v2m[u]] + 1 <=MaxCnt)
		// {
			CNT[v2m[u]]++;
			kclique(l - 1, g, n); // Gv = g[l-1]
			CNT[v2m[u]]--;
		// }

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
	return;
}

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

int main(int argc, char** argv) {
	cout << getTime() << endl;
	auto mp = get_index_mem();
	auto m1 = mp["pk"];

	specialsparse* g;
	partiteSize = atoi(argv[1]);
	unsigned char k = atoi(argv[2]);
	// partiteSize = 4;
	// unsigned char k = 5;
	K = k;

	cout << "# K = " << K << endl;
	unsigned long long n;

	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;
	MaxCnt = K - partiteSize + 1;
	// char dataset[100] = "../../Dataset/emailZJ4.csv";
	// // char dataset[100] = "../data/test";
	// g = readedgelist(dataset);
	printf("Reading edgelist from file %s\n", argv[3]);
	g = readedgelist(argv[3]);
	CNT = new int[partiteSize];
	memset(CNT, 0, partiteSize * sizeof(int));
	// A = atoi(argv[4]);
	// B = atoi(argv[5]);
	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;
	mp = get_index_mem();
	auto graphSize = mp["pk"];
	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);
	printf("Number of partiteCount = %u\n", partiteCount.size());

	printf("Building the graph structure\n");

	mkspecial(g, k);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;
	auto ti= Get_Time();
	printf("Iterate over all cliques\n");

	//qsort_r(g->sub[k],g->ns[k],sizeof(unsigned),compare_r,g->d[k]);
	//arg_bucket_sort(g->sub[k], g->ns[k], g->d[k]);
	n = 0;

	visCNT.resize(g->n);
	C.resize(g->n);
	st.resize(g->n);
	kclique(k, g, &n);

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	mp = get_index_mem();
	auto m2 = mp["pk"];
	fstream fs = fstream("../base_res.csv", ios::app);
	fs << argv[3] << ",base," << g->n << "," << g->e << "," << (int)k << "," << n << "," << Duration(ti) << "," << m2 - m1 << "," << graphSize - m1 << "," << A << "," << B << endl;

	freespecialsparse(g, k);

	printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}