#include <bits/stdc++.h>
using namespace std;

const int T = 1e7 + 5e6;

struct node {
    int l,r,val;
    char lazy;
};

node bag[T];
int nodes;

int a,b,v,n,q;

void push(node *at, int val, int i, int j) {
    at->lazy = val;
    at->val = (at->lazy == 1? (j-i+1) : 0);
}

void prop(node *at, int i, int j, int mid) {
    if(!at->lazy) return;
    push(bag+at->l,at->lazy,i,mid);
    push(bag+at->r,at->lazy,mid+1,j);
    at->lazy = 0;
}

void update(node *at, int i, int j) {
    if(i >= a and j <= b) push(at,v,i,j);
    else {
        int mid = (i+j) >> 1;
        if(!at->l) at->l = ++nodes, at->r = ++nodes;
        prop(at,i,j,mid);
        if(mid >= a) update(bag+at->l,i,mid);
        if(mid < b) update(bag+at->r,mid+1,j);
        at->val = bag[at->l].val + bag[at->r].val;
    }
}

int main() {
    scanf("%d %d", &n, &q);
    node *root = bag;
    while(q--) {
        scanf("%d %d %d", &a, &b, &v);
        update(root,1,n);
        printf("%d ", n-root->val);
    }
    fflush(stdout);
    return 0;
}

