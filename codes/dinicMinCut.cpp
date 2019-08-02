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
 
