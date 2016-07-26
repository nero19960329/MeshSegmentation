#pragma once

class DisjointSet {
private:
    int *p;
    int *rank;
    int size;

public:
    DisjointSet() : p(NULL), rank(NULL), size(0) {}
    DisjointSet(int len) : size(len) {
        p = new int[size];
        rank = new int[size];
        memset(p, 0, size * sizeof(int));
        memset(rank, 0, size * sizeof(int));
    }
    ~DisjointSet() {
        if (p) {
            delete[] p;
            delete[] rank;
        }
    }

    void MakeSet(int x) {
        p[x] = x;
        rank[x] = 0;
    }

    void Union(int x, int y) {
        link(FindSet(x), FindSet(y));
    }

    int FindSet(int x) const {
        if (x != p[x]) {
            p[x] = FindSet(p[x]);
        }
        return p[x];
    }

private:
    void link(int x, int y) {
        if (rank[x] > rank[y]) {
            p[y] = x;
        } else {
            p[x] = y;
            if (rank[x] == rank[y]) {
                ++rank[y];
            }
        }
    }
};