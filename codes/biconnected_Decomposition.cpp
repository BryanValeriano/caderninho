#include <bits/stdc++.h>
using namespace std;

#define pb push_back
#define mk make_pair
#define fi first
#define se second
#define eb emplace_back

typedef long long ll;
typedef pair<int,int> ii;
typedef vector< pair<int,int> > vii;
const int INF = 0x3f3f3f3f;

const int T = 5e4 + 10;
const int N = 1e5;

struct edge { 
    int i,j,id;

    int dif(int x) {
        if(i != x) return i;
        return j;
    }
};

vector<int> g[T];
edge e[T];
vector<int> comps[T];
int low[T];
int tin[T];
bool isCut[T];
bool vis[T];
int t;
int stk[T];

void init() {
    for(int i = 0; i < T; i++) {
        isCut[i] = vis[i] = false;
        tin[i] = low[i] = 0;
        g[i].clear(); 
        comps[i].clear();
    }
    t = 0;
}

void makeComp(int at) {
    set<int> cuts; 
    set<int> tam;
    while(stk[stk[0]] != at) {
        int u = e[stk[stk[0]]].i;
        int v = e[stk[stk[0]]].j;
        if(!isCut[u]) tam.insert(u);
        else cuts.insert(u);
        if(!isCut[v]) tam.insert(v);
        else cuts.insert(v);
        --stk[0];
    }
    if(at != 0) {
        int u = e[stk[stk[0]]].i;
        int v = e[stk[stk[0]]].j;
        if(!isCut[u]) tam.insert(u);
        else cuts.insert(u);
        if(!isCut[v]) tam.insert(v);
        else cuts.insert(v);
        --stk[0];
    }
    if(cuts.size() == 1) comps[*cuts.begin()].pb(tam.size());
}

int dfs(int u, int p, int id) {

    if(!vis[id]) vis[id] = 1, stk[++stk[0]] = id;

    //back edge
    if(tin[u] != 0) {
        low[p] = min(low[p], tin[u]);
        return low[p];
    }

    tin[u] = low[u] = ++t;
    bool hasC = false;
    for(int x : g[u]) {
        if(e[x].id == id) continue;
        
        if(dfs(e[x].dif(u), u, e[x].id) < 0) { 

            low[u] = min(low[u], low[e[x].dif(u)]);

            if(u != p ? low[e[x].dif(u)] >= tin[u] : hasC) { 
                isCut[u] = true; 
                makeComp(e[x].id);
            }
        }

        hasC = true;
    }

    //forth edge
    return -1;
}

ll bin[N][3];

ll choose() { 
   for(int i = 0; i < N; i++) bin[i][0] = 1; 
   for(int j = 0; j <= 2; j++) bin[j][j] = 1;

   for(int i = 1; i < N; i++)
       for(int j = 1; j < i and j <= 2; j++)
            bin[i][j] = bin[i-1][j-1] + bin[i-1][j];
}

int main() {
    ios::sync_with_stdio(false);
    int m, u, v;
    cin >> m;
    choose();
    int z = 0;
    while(m) {
        int cont = 0;
        ll n = 0;
        init();
        for(int i = 0; i < m; i++) {
            cin >> u >> v;
            e[++cont] = {u,v,cont};
            if(g[u].size() == 0) n++;
            g[u].pb(cont);
            if(g[v].size() == 0) n++;
            g[v].pb(cont);
        }
        vis[0] = true;
        dfs(1,1,0);
        if(stk[0] > 0) makeComp(0);
        ll ans = 0;
        ll prod = 1;
        for(int i = 1; i < T; i++) {
            if(comps[i].size() > 0) {
                for(auto x : comps[i]) {
                    ans++;
                    prod *= x;
                }
            }
        }
        if(ans == 0) { ans = 2; prod = bin[n][2]; }  
        cout << "Case " << ++z << ": " << ans << " " << prod << endl;
        stk[0] = 0;
        cin >> m; 
    }
    return 0;
}

