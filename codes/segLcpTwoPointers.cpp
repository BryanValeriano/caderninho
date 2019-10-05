#include <bits/stdc++.h>
using namespace std;

#define fr(i,n) _back
#define pb	push_back
#define eb	emplace_back
#define mk	make_pair
#define fi	first
#define se	second
#define endl	'\n'

typedef long long ll;
typedef pair<ll,ll> pi;
typedef vector<pi> vii;
const ll INF = 0x3f3f3f3f;
const double PI = acos(-1.0);

const int T = 1e5 + 10;
string v[T];
ll ve[T];
ll seg[4*T];
vii pref[4*T];
vii suf[4*T];
ll a,b;
int n;

ll lcp(int x, int y, int w, int z) {
    for(int i = 0; i < min(x,y); i++)
        if(v[w][i] != v[z][i]) return i;
    return min(x,y);
}

void buildCord(int node, int ii, int jj) {
    int esq = 2*node;
    int dir = esq+1;
    pref[node].clear();
    suf[node].clear();

    int start = 0;

    for(pi lef : pref[esq]) pref[node].eb(lef);
    while(start < pref[dir].size() and pref[dir][start].fi >= pref[node].back().fi)
        pref[node].back().se += pref[dir][start++].se;
    while(start < pref[dir].size()) pref[node].eb(pref[dir][start++]);

    start = 0;

    for(pi lef : suf[dir]) suf[node].eb(lef);
    while(start < suf[esq].size() and suf[esq][start].fi >= suf[node].back().fi)
        suf[node].back().se += suf[esq][start++].se;
    while(start < suf[esq].size()) suf[node].eb(suf[esq][start++]);
}

ll twoP(int node, ll i, ll j, ll bordaI, ll bordaJ) {
    ll ans = 0;
    ll x = 0;
    ll qtdL = 0;
    ll qtdR = 0;
    int esq = 2*node;
    int dir = esq+1;
    int mid = (i+j) >> 1;
    bordaI = mid - bordaI + 1;
    bordaJ = bordaJ - mid;

    auto l = suf[esq].begin();
    auto r = pref[dir].begin();
    qtdL += l->se;
    qtdR += r->se;

    while(1) {
        x = min(l->fi,r->fi);
        qtdL = min(bordaI, qtdL);
        qtdR = min(bordaJ, qtdR);

        ans = max(ans,x*(qtdL+qtdR+1));

        auto l1 = l;
        auto r1 = r;
        if(l != suf[esq].end()) l1++;
        if(r != pref[dir].end()) r1++;

        if(qtdR == bordaJ and qtdL == bordaI) break;

        if(qtdR >= bordaJ or (qtdL < bordaI and l1->fi >= r1->fi and l != l1)) qtdL += l1->se, l++;
        else if(qtdR < bordaJ and r != r1) qtdR += r1->se, r++;
        else break;
    }

    return ans;

}

void build(int node, ll i, ll j) {
    if(i == j) {
        seg[node] = max({ve[i]*2,(ll)v[i].size(),(ll)v[i+1].size()});
        pref[node].eb(ve[i],1);
        suf[node].eb(ve[i],1);
    } else {
        int mid = (i+j) >> 1;
        build(2*node,i,mid);
        build(2*node+1,mid+1,j);
        buildCord(node,i,j);
        ll y = twoP(node,i,j,i,j);
        seg[node] = max({y, seg[2*node], seg[2*node+1]});
    }
}

void update(int node, ll i, ll j) {
    if(i > b or j < a) return;
    if(i >= a and j <= b) {
        seg[node] = max({ve[i]*2,(ll)v[i].size(),(ll)v[i+1].size()});
        pref[node] = {mk(ve[i],1)};
        suf[node] = {mk(ve[i],1)};
    } else {
        int mid = (i+j) >> 1;
        update(2*node,i,mid);
        update(2*node+1,mid+1,j);
        buildCord(node,i,j);
        ll y = twoP(node,i,j,i,j);
        seg[node] = max({y, seg[2*node], seg[2*node+1]});
    }
}

ll query(int node, ll i, ll j) {
    if(i > b or j < a) return 0;
    if(i >= a and j <= b) return seg[node];
    int mid = (i+j) >> 1;
    ll w = query(2*node,i,mid);
    ll z = query(2*node+1,mid+1,j);
    ll y = (a <= mid and b > mid? twoP(node,i,j,max(a,i),min(b,j)) : 0);
    return max({w,z,y});
}

int main() {
    ios_base::sync_with_stdio(false);
    int q,op; cin >> n >> q;
    for(int i = 1; i <= n; i++) cin >> v[i];
    for(int i = 1; i < n; i++) ve[i] = lcp(v[i].size(),v[i+1].size(),i,i+1);
    if(n > 1) build(1,1,n-1);

    while(q--) {
        cin >> op >> a;
        if(op == 1) {
            cin >> b;
            if(a == b) cout << v[a].size() << endl;
            else b--, cout << query(1,1,n-1) << endl;
        } else {
            cin >> v[a];
            ve[a] = lcp(v[a].size(), v[a+1].size(), a,a+1);
            b = a;
            if(n > 1) update(1,1,n-1);
            if(a > 1) {
                ve[a-1] = lcp(v[a-1].size(), v[a].size(), a-1,a);
                b--; a--;
                if(n > 1) update(1,1,n-1);
            }
        }
    }

    return 0;
}

