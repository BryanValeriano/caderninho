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

const int T = 1e5 + 3;
vector<int> c(T), p(T), cn(T), pn(T), cnt(T);
string s;
int q;

void sortShifts() {
    int n = s.size();
    int tam = 256;
    fill(cnt.begin(), cnt.begin() + n, 0);

    for(int i = 0; i < n; i++) cnt[s[i]]++;
    for(int i = 1; i < tam; i++) cnt[i] += cnt[i-1];
    for(int i = 0; i < n; i++) p[--cnt[s[i]]] = i;

    int cls = 1;
    c[p[0]] = 0;

    for(int i = 1; i < n; i++) {
        if(s[p[i]] != s[p[i-1]]) cls++;
        c[p[i]] = cls-1;
    }

    int h = 1;

    while(cls != n) {

        for(int i = 0; i < n; i++) {
            pn[i] = p[i] - h;
            if(pn[i] < 0) pn[i] += n;
        }

        fill(cnt.begin(), cnt.begin() + n, 0);

        for(int i = 0; i < n; i++) cnt[c[pn[i]]]++;
        for(int i = 1; i < cls; i++) cnt[i] += cnt[i-1];
        for(int i = n-1; i >= 0; i--) p[--cnt[c[pn[i]]]] = pn[i];

        cn[p[0]] = 0;
        cls = 1;

        for(int i = 1; i < n; i++) {
            ii prev = {c[p[i-1]], p[i-1] + h < n? c[p[i-1] + h] : -1};
            ii at = {c[p[i]], p[i] + h < n? c[p[i] + h] : -1};
            if(at != prev) ++cls;
            cn[p[i]] = cls-1;
        }

        c.swap(cn);
        h <<= 1;
    }
}

void construct() {
    s += "$";
    sortShifts();
    p.erase(p.begin());
    s.pop_back();
}

bool bigger(int pos, string x) {
    for(int i = 0; i < x.size(); i++) {
        if(i+p[pos] >= p.size() or s[p[pos]+i] < x[i]) return false;
        if(s[p[pos]+i] > x[i]) return true;
    }
    return true;
}

bool query(string x) {
    int l = 0;
    int r = s.size()-1;

    while(l <= r) {
        int mid = (l+r) >> 1;
        if(bigger(mid,x)) r = mid-1;
        else l = mid+1;
    }

    for(int i = 0; i < x.size(); i++)
        if(i+p[l] >= p.size() or x[i] != s[p[l]+i]) return false;

    return true;
}

int main() {
    ios::sync_with_stdio(false);

    while(cin >> s) {
        cin >> q;
        construct();
        while(q--) {
            string x;
            cin >> x;
            cout << (query(x)? 'Y' : 'N') << endl;
        }
    }

    return 0;
}

