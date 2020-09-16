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

vector<int> c(T), p(T), cn(T), pn(T), cnt(T);
vector<int> lcp(T), tam(T), rk(T);
string s;

void shift() {
    int n = s.size();
    int alphabet = 256;

    for(int i = 0; i < n; i++) cnt[s[i]]++;
    for(int i = 1; i < alphabet; i++) cnt[i] += cnt[i-1];
    for(int i = 0; i < n; i++) p[--cnt[s[i]]] = i;

    int cls = 1;
    c[p[0]] = 0;

    for(int i = 1; i < n; i++) {
        cls += (s[p[i]] != s[p[i-1]]);
        c[p[i]] = cls-1;
    }

    int h = 1;

    while(h < n) {
        for(int i = 0; i < n; i++) {
            pn[i] = p[i] - h;
            if(pn[i] < 0) pn[i] += n;
        }

        fill(cnt.begin(), cnt.begin() + cls, 0);

        for(int i = 0; i < n; i++) cnt[c[pn[i]]]++;
        for(int i = 1; i < cls; i++) cnt[i] += cnt[i-1];
        for(int i = n-1; i >= 0; i--) p[--cnt[c[pn[i]]]] = pn[i];

        cls = 1;
        cn[p[0]] = 0;

        for(int i = 1; i < n; i++) {
            ii prev = {c[p[i-1]], p[i-1] + h < n? c[p[i-1]+h] : -1};
            ii at = {c[p[i]], p[i] + h < n? c[p[i]+h] : -1};
            cls += (at != prev);
            cn[p[i]] = cls-1;
        }

        h <<= 1;
        //if(cls == n) break;
        c.swap(cn);
    }
}

void longestCp() {
    int n = s.size();
    for(int i = 0; i < n; i++) rk[p[i]] = i;
    int k = 0;
    for(int i = 0; i < n; i++) {
        if(rk[i] == n-1) { k = 0; continue; }
        int j = p[rk[i]+1];
        while(i+k < n && j+k < n && s[i+k] == s[j+k]) k++;
        lcp[rk[i]] = k;
        if(k) k--;
    }
}

void buildTam() {
    int n = s.size();
    tam[0] = n-p[0];
    for(int i = 1; i < n; i++) tam[i] = n-p[i]-lcp[i-1];
    for(int i = 1; i < n; i++) tam[i] += tam[i-1];
}

string search(int x) {
    int n = s.size();
    int y = lower_bound(tam.begin(), tam.begin()+n,x) - tam.begin();

    string ans;
    int i = p[y];

    if(y > 0) {
        x -= tam[y-1];
        while(i-p[y] < lcp[y-1]) ans += s[i++];
    }

    while(x) ans += s[i++], x--;
    ans += "\n";
    return ans;
}


int main() {
    ios_base::sync_with_stdio(false);
    cin >> s;
    buildSuffix();
    longestCp();
    buildTam();

    int q, x; cin >> q;
    string ans;
    while(q--) {
        cin >> x;
        ans += search(x);
    }
    cout << ans;
    return 0;
}

