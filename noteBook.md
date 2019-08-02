# Trick Implementations



### Biconnected Decomposition

```c++
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

```







### moHilbert

```c++
#include <bits/stdc++.h>
using namespace std;

#define pb push_back
#define mk make_pair
#define fi first
#define se second
#define eb emplace_back

typedef long long ll;
typedef unsigned long long ull;
typedef pair<int,int> ii;
typedef vector< pair<int,int> > vii;
const int INF = 0x3f3f3f3f;

const int pw = 21; 
const int T = 1e5 + 1000;
const int N = sqrt(T);

ull hilbert(int x, int y, int pow, int rotate) {
    if(!pow) return 0;
    int hpow = 1 << (pow-1);
    int seg = (x < hpow) ? (
            (y < hpow) ? 0 : 3
        ) : (
            (y < hpow) ? 1 : 2
        );
    seg = (seg + rotate) & 3;
    const int rotateDelta[4] = {3,0,0,1};
    int nx = x & (x ^ hpow), ny = y & (y ^ hpow);
    int nroot = (rotate + rotateDelta[seg]) & 3;
    ull subSquareSize = ull(1) << (2*pow - 2);
    ull ans = seg * subSquareSize;
    ull add = hilbert(nx,ny,pow-1, nroot);
    ans += (seg == 1 or seg == 2) ? add : (subSquareSize - add - 1);
    return ans;
}

struct query { 
    int ind, l, r, k;
    ull ord;
    void getO() {
        ord = hilbert(l,r,pw,0);
    }
    bool operator < (const query &b) const {
        return ord < b.ord;
    }
};

vector<query> qrs;
vector<int> g[T];
int cor[T];
int tmp[T];
int tin[T];
int tout[T];
int resp[T];
int track[T];
int ans[T];
bool vis[T];
int n, q, t;

void add(int ind) {
    track[cor[ind]]++;
    resp[track[cor[ind]]]++;
}

void tira(int ind) {
    resp[track[cor[ind]]]--;
    track[cor[ind]]--;
}

void mo() {
    int l = 1;
    int r = 1;
    add(1);
    for(auto x : qrs) {
        while(r < x.r) add(++r);
        while(r > x.r) tira(r--);
        while(l < x.l) tira(l++);
        while(l > x.l) add(--l);
        ans[x.ind] = resp[x.k];
    }
}

void dfs(int u) {
    vis[u] = true;
    tin[u] = ++t;
    cor[t] = tmp[u];
    for(int v : g[u]) 
        if(!vis[v]) dfs(v);
    tout[u] = t;
}

int main() {
    ios::sync_with_stdio(false);
    cin >> n >> q;
    int u,v;
    for(int i = 1; i <= n; i++) cin >> tmp[i];
    for(int i = 0; i < n-1; i++) {
        cin >> u >> v;
        g[u].pb(v);
        g[v].pb(u);
    }
    dfs(1);
    for(int i = 0; i < q; i++) {
        cin >> u >> v;
        qrs.pb({i,tin[u], tout[u],v,0});
        qrs.rbegin()->getO();
    }
    sort(qrs.begin(), qrs.end());
    mo();
    for(int i = 0; i < q; i++) cout << ans[i] << endl;
    return 0;
}

```





### Dinic MinCut

```c++
#include <bits/stdc++.h>
using namespace std;
 
#define pb push_back
#define eb emplace_back
#define mk make_pair
#define fi first
#define se second
#define cc(x)	cout << #x << " = " << x << endl
#define ok		cout << "ok" << endl
#define endl '\n'
 
typedef long long ll;
typedef pair<int,int> ii;
const int INF = 0x3f3f3f3f;
const double PI = acos(-1.0);
 
 
const int N = 210;
string id[N];

struct edge {
	int a, b, cap, flow;
};
 
int d[N], ptr[N], q[N];
vector<edge> e;
vector<int> g[N];
 
void add(int a, int b, int cap) {
	edge e1 = { a, b, cap, 0 };
	edge e2 = { b, a, 0, 0 };
	g[a].push_back ((int) e.size());
	e.push_back (e1);
	g[b].push_back ((int) e.size());
	e.push_back (e2);
}
 
bool bfs(int s, int t) {
	int qh=0, qt=0;
	q[qt++] = s;
	memset (d, -1, sizeof d);
	d[s] = 0;
	while (qh < qt && d[t] == -1) {
		int v = q[qh++];
		for (size_t i=0; i<g[v].size(); ++i) {
			int id = g[v][i],
				to = e[id].b;
			if (d[to] == -1 && e[id].flow < e[id].cap) {
				q[qt++] = to;
				d[to] = d[v] + 1;
			}
		}
	}
	return d[t] != -1;
}
 
int dfs (int v, int flow, int t) {
	if (!flow)  return 0;
	if (v == t)  return flow;
	for (; ptr[v]<(int)g[v].size(); ++ptr[v]) {
		int id = g[v][ptr[v]],
			to = e[id].b;
		if (d[to] != d[v] + 1)  continue;
		int pushed = dfs (to, min (flow, e[id].cap - e[id].flow), t);
		if (pushed) {
			e[id].flow += pushed;
			e[id^1].flow -= pushed;
			return pushed;
		}
	}
	return 0;
}
 
int dinic(int s, int t) {
	int flow = 0;
    while(bfs(s, t)) { 
		memset (ptr, 0, sizeof ptr);
		while (int pushed = dfs (s, INF, t)) flow += pushed;
	}
	return flow;
}
 
ii dig(string x) {
	ii ans = mk(0,0);
	for(int i = 0; i < 3; i++) ans.fi += (x[i] - '0');
	for(int i = 3; i < 6; i++) ans.se += (x[i] - '0');
	return ans;
}
 
int vis[N];
 
void find_MinCut(int at, int s) {
	vis[at] = 1;
    if(at >= 102) return;
 
	for(int i = 0; i < g[at].size(); i++) {
		edge edg = e[g[at][i]];
		int next = edg.b;
		if(edg.flow == edg.cap and next != s and at != s and next >= 102 and next < N-1) {
			ii a = dig(id[at]);
			ii b = dig(id[next]);
			if(a.fi == b.se) cout << "AT " << id[at] << " " << id[next] << endl;
			else cout << "TA " << id[next] << " " << id[at] << endl;
		} 
        else if(!vis[next]) find_MinCut(next, s);
	}
}
 
int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
 
	ios_base::sync_with_stdio(false);
	int n, m;
	cin >> n >> m;
 
	for(int i = 1; i <= n; i++) {
		cin >> id[i];
		add(0,i,1);
	}
	for(int i = 102; i < 102+m; i++) {
		cin >> id[i]; 
		add(i, N-1, 1);
	}
 
	for(int i = 1; i <= n; i++) {
		for(int j = 102; j < 102+m; j++) {
			ii a = dig(id[i]);
			ii b = dig(id[j]);
			if(a.fi == b.se or a.se == b.fi) add(i,j,1);
		}
	}
	
	cout << dinic(0, N-1) << endl;
	find_MinCut(0,0);
 
	return 0;
}
```





### Mat exp

```c++
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
const int MOD = 1e9 + 9;
const int n = 3;

typedef vector< vector<ll> > mat;

mat operator * (const mat &a, const mat &b) {
    mat m(n, vector<ll>(n));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) 
            for(int k = 0; k < n; k++) {
                m[i][j] += (a[i][k] * b[k][j]) % MOD;
                m[i][j] %= MOD;
            }
    return m;
}

mat id() {
    mat m(n, vector<ll>(n));
    for(int i = 0; i < n; i++) m[i][i] = 1;
    return m;
}

mat trans() {
    mat m(n, vector<ll>(n,1));
    m[1][1] = m[1][2] = m[2][0] = m[2][2] = 0;
    return m;
}

void print(mat m) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++)
            cout << m[i][j] << " ";
        cout << endl;
    }
}

ll expo(ll e) {
    mat ans = id();
    mat ele = trans();

    while(e) {
        if(e & 1) ans = ans * ele;
        ele = ele * ele;
        e >>= 1;
    }
    
    ll resp;
    resp = (2 * ans[0][0]) % MOD;
    resp = (resp + ans[0][1]) % MOD;
  
    return resp;
}


int main() {
    ios::sync_with_stdio(false);
    ll p;
    cin >> p;
    while(p != 0) {
        if(p == 1) cout << 0 << endl;
        else if(p == 2) cout << 1 << endl;
        else if(p == 3) cout << 2 << endl;
        else cout << expo(p-3) << endl;
        cin >> p;
    }
    return 0;
}

```

