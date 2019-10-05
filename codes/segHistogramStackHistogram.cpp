#include <bits/stdc++.h>
using namespace std;

#define pb push_back
#define mk make_pair
#define fi first
#define se second
#define eb emplace_back
#define endl '\n'

typedef long long ll;
typedef pair<int,int> ii;
typedef vector< pair<int,int> > vii;
const int INF = 0x3f3f3f3f;

const int N = 1e5 + 6;
string v[N];
int ve[N];
const int T = 78;
int seg[T + 3][4*N + 3];
int suf[T + 3][4*N + 3];
int pref[T + 3][4*N + 3];
int corn[4*N+3];
int n, a, b, val;
set<int> big;

int lcp(int x, int y) {
    for(int i = 0; i < min(v[x].size(),v[y].size()); i++)
        if(v[x][i] != v[y][i]) return i;
    return min(v[x].size(), v[y].size());
}

void upd(int node, int i, int j) {
    if(i > b or j < a) return;
    if(i >= a and j <= b) {
        corn[node] = v[i].size();
        return;
    }
    int mid = (i+j) >> 1;
    upd(2*node,i,mid);
    upd(2*node+1,mid+1,j);
    corn[node] = max(corn[2*node],corn[2*node+1]);
}

int qer(int node, int i, int j) {
    if(i > b or j < a) return 0;
    if(i >= a and j <= b) return corn[node];
    int mid = (i+j) >> 1;
    return max(qer(2*node,i,mid),qer(2*node+1,mid+1,j));
}

void update(int h, int node, int i, int j) {
    if(i > b or j < a) return;
    if(i >= a and j <= b) {
        seg[h][node] += val;
        pref[h][node] += val;
        suf[h][node] += val;
    } else {
        int mid = (i+j) >> 1;
        int l = 2*node;
        int r = l+1;
        update(h, l, i, mid);
        update(h, r, mid+1, j);
        suf[h][node] = suf[h][r] + (suf[h][r] == (j-mid)? suf[h][l] : 0);
        pref[h][node] = pref[h][l] + (pref[h][l] == (mid-i+1)? pref[h][r] : 0);
        seg[h][node] = max({seg[h][r],seg[h][l],suf[h][l]+pref[h][r]});
    }
}

ll query(int h, int node, int i, int j) {
    if(i > b or j < a) return 0;
    if(i >= a and j <= b) return seg[h][node];
    int mid = (i+j) >> 1;
    int l = 2*node;
    int r = l+1;
    ll y = (a <= mid and b > mid)? min(max(0,mid-a+1), suf[h][l])  + min(max(0,b-mid), pref[h][r]) : 0;
    return max({query(h,l,i,mid), query(h,r,mid+1,j), y});
}

ll stackWay() {
    vector<ll> bag;
    int last = -1;

    for(auto x : big) {
        if(x >= a and x <= b) {
            if(last != x-1) bag.pb(0);
            bag.pb(ve[x]);
            last = x;
        }
    }

    ll i = 0;
    stack<ll> s;
    ll ans = 0;

    while(i < bag.size()) {
        if(s.empty() or bag[i] >= bag[s.top()]) s.push(i++);
        else {
            int st = s.top();
            s.pop();
            ans = max(ans, bag[st] * ((s.empty()? i : i-s.top()-1) + 1));
        }
    }

    while(!s.empty()) {
        int st = s.top();
        s.pop();
        ans = max(ans, bag[st] * ((s.empty()? i : i-s.top()-1) + 1));
    }

    return ans;
}

void upCol(int id, int z) {
    val = z;
    for(int j = 1; j <= T; j++) {
        if(ve[id] >= j) a = b = id, update(j,1,1,n);
        else break;
    }

    if(ve[id] > T) {
        if(z == 1) big.insert(id);
        else big.erase(id);
    }
}

ll subAlg() {
    ll ans = 0;
    for(ll j = 1; j <= T; j++) {
        ll w = query(j,1,1,n);
        if(w) ans = max(ans,j*(w+1));
        else break;
    }
    ans = max(ans,stackWay()); 
    b++;
    ans = max(ans,(ll)qer(1,1,n));
    return ans;
}
 
int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);
    int q;
    cin >> n >> q;
    for(int i = 1; i <= n; i++) cin >> v[i];
    for(int i = 1; i < n; i++) ve[i] = lcp(i,i+1);
    for(int i = 1; i < n; i++) upCol(i,1);
    for(int i = 1; i <= n; i++) a = b = i, upd(1,1,n);

    int op;
    int x,y;

    while(q--) {
        cin >> op >> x;
        if(op == 1) {
            cin >> y;
            if(x == y) cout << v[x].size() << endl;
            else {
                a = x;
                b = y-1;
                cout << subAlg() << endl;
            }
        } else {
            cin >> v[x];
            a = b = x;
            upd(1,1,n);
            upCol(x,-1);
            ve[x] = lcp(x,x+1);
            upCol(x,1);
            if(x > 1) {
                upCol(x-1,-1);
                ve[x-1] = lcp(x-1,x);
                upCol(x-1,1);
            }
        }
    }
    
    return 0;
}

