#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>

#include <mst.h>

typedef long long lld;
typedef unsigned long long llu;
typedef unsigned int uint;
using namespace std;

struct pq_entry
{
    double x;
    pair<int, int> e;
    bool operator <(const pq_entry &a) const
    {
        return (x > a.x);
    }
};

inline vector<int> prim(int n, double adj[MAX_N][MAX_N])
{
    vector<bool> mark(n + 1, false);
    vector<int> ret(n + 1, -1);
    priority_queue<pq_entry> pq_prim;
    int xt = 1;
    int amt = 0;
    while (amt < n - 1)
    {
        mark[xt] = true;
        for (int i=1;i<=n;i++)
        {
            if (!mark[i]) pq_prim.push({adj[max(xt, i)][min(xt, i)], {xt, i}});
        }
        while (!pq_prim.empty())
        {
            auto X = pq_prim.top();
            pq_prim.pop();
            if (!mark[X.e.second])
            {
                ret[X.e.second] = X.e.first;
                xt = X.e.second;
                amt++;
                break;
            }
        }
    }
    return ret;
}

inline pair<int, vector<vector<int> > > build_tree(vector<int> p)
{
    int root;
    vector<vector<int> > T(p.size(), vector<int>());
    for (uint i=1;i<p.size();i++)
    {
        if (p[i] == -1) root = i;
        else T[p[i]].push_back(i);
    }
    return {root, T};
}

void preorder(int x, vector<vector<int> > T, vector<int> &ret)
{
    ret.push_back(x);
    for (uint i=0;i<T[x].size();i++)
    {
        preorder(T[x][i], T, ret);
    }
}

pair<double, vector<pair<int, int> > > tsp_mst(int n, double adj[MAX_N][MAX_N])
{
    vector<int> mst = prim(n, adj);
    auto mst_t = build_tree(mst);
    vector<int> walk;
    preorder(mst_t.first, mst_t.second, walk);
    
    double val = 0.0;
    vector<pair<int, int> > used;
    
    for (uint i=0;i<walk.size();i++)
    {
        int x = walk[i];
        int y = walk[(i + 1) % walk.size()];
        val += adj[max(x, y)][min(x, y)];
        used.push_back({x, y});
    }
    
    return {val, used};
}
