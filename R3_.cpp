//
// Created by gaga on 25-2-24.
//
#include <bits/stdc++.h>
using namespace std;
#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast< \
std::chrono::microseconds>(Get_Time() - start).count() / (float)1000
#define Print_Time(str, start) std::cout << str << Duration(start) << " ms" << std::endl
#define bzero(a, b)             memset(a, 0, b)
#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed
#define DEBUG 0
using namespace std;
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
bool *res;
bool *vis;
bool *set_vis;
unsigned *ans;
unsigned f;
void readFile(string fileName)
{
    std::ifstream input(fileName);
    if (!input) {
        std::cerr << "无法打开文件: " << fileName << std::endl;
        return;
    }
    input >> g->n >> g->e;
    //cout << "g->n = " << g->n << endl;
    res = new bool[g->n];
    belong = new unsigned[g->n];
    vis = new bool[g->n];
    g->d = new unsigned[g->n];
    g->edges = new edge[g->e];
    g->cd = new unsigned[g->n+1];
    g->adj = new unsigned[2 * g->e];
    for (int i=0; i<g->n; i++)
    {
        input >> belong[i];
        f = max(f, belong[i]);
    }
    f ++;
    set_vis = new bool[f];
    ans = new unsigned[f];
    for (int i=0; i<g->e; i++)
    {
        input >> g->edges[i].s >> g->edges[i].t;
        g->d[g->edges[i].s] ++;
        g->d[g->edges[i].t] ++;
    }
    g->cd[0] = 0;
    for (int i=1; i<=g->n; i++)
    {
        g->cd[i] = g->cd[i-1] + g->d[i-1];
        g->d[i-1] = 0;
    }
    for (int i=0; i<g->e; i++)
    {
        unsigned u = g->edges[i].s, v = g->edges[i].t;
        g->adj[g->cd[u] + g->d[u]++] = v;
        g->adj[g->cd[v] + g->d[v]++] = u;
    }
    for (int i=0; i<g->n; i++)
    {
        sort(g->adj+g->cd[i], g->adj+g->cd[i+1]);
    }
#if  DEBUG
    for (int i=0; i<g->n; i++)
    {
        cout << i << ":";
        for (int j=g->cd[i]; j<g->cd[i+1]; j++)
        {
            cout << g->adj[j] << " ";
        }
        cout << "\n";
    }
#endif
    input.close();
}
bool ok;
void dfs(unsigned now, unsigned cnt)
{
    if (ok) return;
    ans[cnt] = now;
    if (cnt == f - 1)
    {
        unsigned *p = lower_bound(g->adj+g->cd[now], g->adj+g->cd[now+1], ans[0]);
        if (p == g->adj+g->cd[now+1]) return;
        // int fff = 0;
        // for (int i=g->cd[now]; i<g->cd[now+1]; i++)
        // {
        //     if (g->adj[i] == ans[0])
        //     {
        //         fff = 1;
        //     }
        // }
        if (*p == ans[0])
        {
            ok = true;
            for (int i=0; i<=cnt; i++)
            {
                res[ans[i]] = true;
            }
#if DEBUG
            for (int i=0; i<=cnt; i++)
            {
                cout << belong[ans[i]] << " ";
            }
            cout << "\n";
#endif
        }

        return;
    }
    if (ok) return;
    vis[now] = true;
    set_vis[belong[now]] = true;
    for (int i=g->cd[now]; i<g->cd[now+1]; i++)
    {
        if (ok) break;
        int to = g->adj[i];
        if (vis[to] || set_vis[belong[to]]) continue;
        dfs(to, cnt+1);
    }
    vis[now] = false;
    set_vis[belong[now]] = false;
}

int main(int argc, char** argv)
{
    g = new specialsparse;
    char dataset[100] = "../../Dataset/youtubeXL6.csv";
    readFile(dataset);

    for (int i=0; i<g->n; i++)
    {
        ok = 0;
        if (res[i]) continue;
        // if (i % 10000 == 0)
        cout << i << endl;
        dfs(i, 0);
    }
    string R2 = dataset;
    R2 += "R3";
    fstream output (R2, fstream::out);
    for (int i=0; i<g->n; i++)
    {
        output << res[i];
        if (i != g->n-1) output << "\n";
    }
}
