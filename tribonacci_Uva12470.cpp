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

