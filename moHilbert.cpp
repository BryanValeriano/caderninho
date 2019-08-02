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

