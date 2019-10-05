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
typedef tuple<int,int,int,int> ti;
typedef vector<ti> vti;
const int INF = 0x3f3f3f3f;

const int T = 1e5 + 10;

int n;
vti v;
int X1, Y1, X2, Y2;
int a,b;
ll val;
ll p;
ll ok;
ll len[4*T];
int cont[4*T];

void update(int node, int i, int j) {
    if(i > b or j < a) return;
    if(i >= a and j <= b) cont[node] += val;
    else if(i != j) {
        int mid = (i+j) >> 1; 
        update(2*node,i,mid);
        update(2*node+1,mid+1,j);
    }
    if(cont[node]) len[node] = (j-i+1);
    else if(i != j) len[node] = (len[2*node] + len[2*node+1]);
    else len[node] = 0;
}

bool check(int r) {
    vti tmp;

    for(int i = 0; i < (int)v.size(); i++) {
        int x1,x2,y1,y2;
        tie(x1,y1,x2,y2) = v[i];
        x1 = max(x1-r,X1), x2 = min(x2+r, X2);
        y1 = max(y1-r+1,Y1+1), y2 = min(y2+r,Y2);
        tmp.eb(x1, y1, y2, 1);
        tmp.eb(x2, y1, y2, -1);
    }

    sort(tmp.begin(), tmp.end());

    ll nok = 0;
    ll last = 0;

    for(auto ys : tmp) {
        ll x;
        tie(x,a,b,val) = ys;
        nok += (x - last) * len[1];
        update(1,Y1,Y2);
        last = x;
    }

    return nok*100ll >= ok*p;
}
 
int main() {
    ios_base::sync_with_stdio(false);
    cin >> n;

    int a,b,c,d;
    for(int i = 0; i < n; i++) {
        cin >> a >> b >> c >> d;
        v.eb(a,b,c,d);
    }

    cin >> p >> X1 >> Y1 >> X2 >> Y2; 
    
    ok = (ll)(X2 - X1)*(ll)(Y2 - Y1);
    int l = 0;
    int r = 1e5 + 10;
    int ans = r;

    while(l <= r) {
        int mid = (l+r) >> 1;
        if(check(mid)) ans = mid, r = mid-1;
        else l = mid+1;
    }

    cout << ans << endl;
    //cout << check(5) << endl;
    return 0;
}

