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

    unsigned n;//number of nodes
    unsigned e;//number of edges
    edge *edges;//list of edges

    unsigned *d;//d[l]: degrees of G_l
    unsigned *cd;//cumulative degree: (start with 0) length=n+1
    unsigned *adj;//list of neighbors with lower degree

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
    if (g->n >= 1) g->d[g->n-1] = 0;
    for (int i=0; i<g->e; i++)
    {
        unsigned u = g->edges[i].s, v = g->edges[i].t;
        g->adj[g->cd[u] + g->d[u]++] = v;
        g->adj[g->cd[v] + g->d[v]++] = u;
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
bool dfs(unsigned now, unsigned cnt)
{
    ans[cnt] = now;
    if (cnt == f - 1)
    {
        unsigned *p = lower_bound(g->adj+g->cd[now], g->adj+g->cd[now+1], ans[0]);
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
                cout << ans[i] << " ";
            }
            cout << "\n";
#endif
            return true;
        }

        return false;
    }
    if (ok) return true;
    vis[now] = true;
    set_vis[belong[now]] = true;
    for (int i=g->cd[now]; i<g->cd[now+1]; i++)
    {
        int to = g->adj[i];
        if (vis[to] || set_vis[belong[to]]) continue;
        if (dfs(to, cnt+1))
            return true;
    }
    vis[now] = false;
    set_vis[belong[now]] = false;

    return false;
}
int main(int argc, char** argv)
{
    g = new specialsparse;
    // readFile("test.txt");
    readFile(argv[1]);
    auto t = Get_Time();
    for (int i=0; i<g->n; i++)
    {
        ok = 0;
        if (res[i]) continue;
        dfs(i, 0);
    }
    string s = argv[1];
    s += "R3";
    ofstream output(s);
    for (int i=0; i<g->n; i++)
    {
        output << res[i];
        if (i != g->n-1) output << "\n";
    }
    Print_Time("All Time: ", t);
#if DEBUG
    for (int i=0; i<g->n; i++)
    {
        cout << "i = " << i << " " << res[i] << "\n";
    }
#endif
}
