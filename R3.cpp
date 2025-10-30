#include <bits/stdc++.h>
// #include <tbb/tbb.h>
using namespace std;
#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast< \
std::chrono::microseconds>(Get_Time() - start).count() / (float)1000
#define Print_Time(str, start) std::cout << str << Duration(start) << " ms" << std::endl
#define NLINKS 100000000
#define DEBUG 0

typedef struct {
    unsigned s;
    unsigned t;
} edge;

typedef struct {
    unsigned n;
    unsigned e;
    edge *edges;
    unsigned *d;
    unsigned *cd;
    unsigned *adj;
} specialsparse;

specialsparse *g;
unsigned *belong;
unsigned f;
bool* pass;
void readFile(string fileName) {
    ifstream input(fileName);
    if (!input) {
        cerr << "无法打开文件: " << fileName << endl;
        return;
    }
    input >> g->n >> g->e;
    belong = new unsigned[g->n];
    g->d = new unsigned[g->n]();
    g->edges = new edge[g->e];
    g->cd = new unsigned[g->n + 1];
    g->adj = new unsigned[2 * g->e];

    for (int i = 0; i < g->n; i++) {
        input >> belong[i];
        f = max(f, belong[i]);
    }
    f++;


    for (int i = 0; i < g->e; i++) {
        input >> g->edges[i].s >> g->edges[i].t;
        g->d[g->edges[i].s]++;
        g->d[g->edges[i].t]++;
    }

    g->cd[0] = 0;
    for (int i = 1; i <= g->n; i++) {
        g->cd[i] = g->cd[i - 1] + g->d[i - 1];
        g->d[i - 1] = 0;
    }

    for (int i = 0; i < g->e; i++) {
        unsigned u = g->edges[i].s, v = g->edges[i].t;
        g->adj[g->cd[u] + g->d[u]++] = v;
        g->adj[g->cd[v] + g->d[v]++] = u;
    }

    for (int i = 0; i < g->n; i++) {
        sort(g->adj + g->cd[i], g->adj + g->cd[i + 1]);
    }
    input.close();
    pass = new bool[g->n];
    memset(pass, false, sizeof(bool) * g->n);
    input = ifstream(fileName+"R2");
    if (!input) {
        cerr << "无法打开文件: R2" << fileName << endl;
        return;
    }
    for (int i = 0; i < g->n; i++)
    {
        input >> pass[i];
    }
}

class Work
{

};

void dfs(unsigned now, unsigned cnt, bool* vis, bool* set_vis, unsigned* ans, bool* res, bool& ok, int* visCNT) {
    if (ok) return;
    ans[cnt] = now;
    if (cnt == f - 1) {
        // unsigned *p = lower_bound(g->adj + g->cd[now], g->adj + g->cd[now + 1], ans[0]);
        // if (p == g->adj + g->cd[now + 1]) return;
        // if (*p == now) {
        if (visCNT[now] == 1) {
            ok = true;
            for (int i = 0; i <= cnt; i++) {
                res[ans[i]] = true;
            }
        }
        return;
    }
    if (ok) return;
    vis[now] = true;
    set_vis[belong[now]] = true;
    for (int i = g->cd[now]; i < g->cd[now + 1]; i++) {
        if (ok) break;
        int to = g->adj[i];
        if (!pass[to]) continue;
        if (vis[to] || set_vis[belong[to]]) continue;
        dfs(to, cnt + 1, vis, set_vis, ans, res, ok, visCNT);
    }
    vis[now] = false;
    set_vis[belong[now]] = false;
}

int main(int argc, char** argv) {
    g = new specialsparse;
    // string dataset = "../Dataset/youtubeXL6.csv";
    // congressZJ2 CTZJ2 SDZJ2 youtubeZJ2
    // congressZJ4 CTZJ4 SDZJ4 youtubeZJ4
    // congressXL5 CTXL5 SDXL5 youtubeXL5
    // congressXL6 CTXL6 SDXL6 youtubeXL6
    char dataset[200] = "../../Dataset/KPC_Dataset/actor.txt";
    // string dataset = "../../Dataset/youtubeZJ4.csv";
    readFile(dataset);

    bool* res = new bool[g->n];
    memset(res, false, sizeof(bool) * g->n);
    auto ti = Get_Time();
    // tbb::parallel_for(tbb::blocked_range<int>(0, g->n), [&](const tbb::blocked_range<int>& r) {
        bool* vis = new bool[g->n];
        bool* set_vis = new bool[f];
        auto* ans = new unsigned[f];
        int* visCNT = new int[g->n];
        memset(vis, false, sizeof(bool) * g->n);
        memset(set_vis, false, sizeof(bool) * f);
        memset(ans, 0, sizeof(unsigned) * f);
        memset(visCNT, 0, sizeof(int) * g->n);
        // for (auto i = r.begin(); i < r.end(); i++) {
        for (int i = 0; i < g->n; i++) {
            bool ok = false;
            if (!pass[i]) continue;
            if (res[i]) continue;
            for (int j = g->cd[i]; j < g->cd[i + 1]; j++)
            {
                int to = g->adj[j];
                visCNT[to]++;
            }
            dfs(i, 0, vis, set_vis, ans, res, ok, visCNT);
            for (int j = g->cd[i]; j < g->cd[i + 1]; j++)
            {
                int to = g->adj[j];
                visCNT[to]--;
            }
        }
        delete[] vis;
        delete[] set_vis;
        delete[] ans;
        delete[] visCNT;
    // }, tbb::static_partitioner());
    Print_Time("", ti);
    string R2 = dataset;
    R2 += "R3";
    ofstream output(R2);
    for (int i = 0; i < g->n; i++) {
        output << 1 << "\n";
        // output << res[i] << "\n";
    }
}
