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
const ll INF = 0x3f3f3f3f3f3f3fll;

const int T = 405;
ll mini[2][T][T][T];
int id[T];
vector<int> down;
vector<int> up;
int n,r;
vii v;

void floyd() {
    for(int i = 1; i <= n; i++)
        for(int j = 0; j <= 1; j++)
            mini[j][0][i][i] = 0;

    for(int k = 1; k <= n; k++)
        for(int i = 1; i <= n; i++)
            for(int j = 1; j <= n; j++) {
                mini[0][0][i][j] = min(mini[0][0][i][k] + mini[0][0][k][j], mini[0][0][i][j]);
                mini[0][k][i][j] = min(mini[0][k][i][j], mini[0][0][i][j]);
            }

    for(int k = n; k >= 1; k--)
        for(int i = 1; i <= n; i++)
            for(int j = 1; j <= n; j++) {
                mini[1][0][i][j] = min(mini[1][0][i][k] + mini[1][0][k][j], mini[1][0][i][j]);
                mini[1][k][i][j] = min(mini[1][k][i][j], mini[1][0][i][j]);
            }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);
    cin >> n >> r;

    for(int k = 0; k <= n+1; k++)
        for(int i = 0; i <= n+1; i++)
            for(int j = 0; j <= n+1; j++)
                for(int p = 0; p <= 1; p++)
                    mini[p][k][i][j] = INF;

    for(int i = 0; i < n; i++) {
        ll x; cin >> x;
        v.eb(x,i+1);
    }

    sort(v.begin(), v.end());

    ll last= -INF;
    for(int i = 0; i < n; i++) {
        id[v[i].se] = i+1;
        if(v[i].fi != last) up.pb(i+1);
        last = v[i].fi;
    }

    last = -INF;
    for(int i = n-1; i >= 0; i--) {
        if(v[i].fi != last) down.pb(i+1);
        last = v[i].fi;
    }

    sort(down.begin(),down.end());
    sort(up.begin(),up.end());

    for(int i = 0; i < r; i++) {
        ll a,b,d;
        cin >> a >> b >> d;
        mini[0][0][id[a]][id[b]] = d;
        mini[0][0][id[b]][id[a]] = d;
        mini[1][0][id[a]][id[b]] = d;
        mini[1][0][id[b]][id[a]] = d;
    }

    floyd();

    int q; cin >> q;
    while(q--) {
        int a,b,k,t;
        cin >> a >> b >> k >> t;
        if(t) {
            k = up.size()-k;
            k = max(0,min(k, (int)up.size()-1));
            k = up[k];
        } else {
            k = max(0, min(k, (int)down.size())-1);
            k = down[k];
        }
        ll ans = mini[t][k][id[a]][id[b]];
        cout << (ans == INF? -1 : ans) << endl;
    }

    return 0;
}

