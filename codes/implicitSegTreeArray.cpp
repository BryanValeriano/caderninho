#include <bits/stdc++.h>
using namespace std;

const int T = 1e7+5e6;
int l[T], r[T], seg[T];
char lazy[T];
int a,b,v,n,q,nodes;

void push(int node, int val, int i, int j) {
    lazy[node] = val;
    seg[node] = (lazy[node] == 1? (j-i+1) : 0);
}

void prop(int node, int i, int j, int mid) {
    if(!lazy[node]) return;
    push(l[node],lazy[node],i,mid);
    push(r[node],lazy[node],mid+1,j);
    lazy[node] = 0;
}

void update(int node, int i, int j) {
    if(i >= a and j <= b) push(node,v,i,j);
    else {
        int mid = (i+j) >> 1;
        if(!l[node]) l[node] = ++nodes, r[node] = ++nodes;
        prop(node,i,j,mid);
        if(mid >= a) update(l[node],i,mid);
        if(mid < b) update(r[node],mid+1,j);
        seg[node] = seg[l[node]] + seg[r[node]];
    }
}

int main() {
    scanf("%d %d", &n, &q);
    while(q--) {
        scanf("%d %d %d", &a, &b, &v);
        update(0,1,n);
        printf("%d\n", n-seg[0]);
    }
    fflush(stdout);
    return 0;
}
