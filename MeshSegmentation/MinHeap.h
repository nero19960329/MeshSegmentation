#pragma once

#include "List.h"

template <class heapElem, class heapElemComp>
class MinHeap {
private:
    heapElem *a;
    int *pos;
    int size;
    heapElemComp cmpFun;

public:
    MinHeap() : a(NULL), size(0) {}

    MinHeap(heapElem *input, int len) : size(len) {
        a = new heapElem[size];
        for (int i = 0; i < size; ++i) {
            a[i] = input[i];
        }

        pos = new int[size];
        for (int i = 0; i < size; ++i) {
            pos[i] = i;
        }

        for (int i = parent(size); i >= 0; --i) {
            minHeapify(i);
        }
    }

    MinHeap(List<heapElem> list, int maxIndex) {
        size = list.Size();

        // still need a map ...

        a = new heapElem[maxIndex];
        pos = new int[maxIndex];
        int i = 0;
        for (List<heapElem>::iterator it = list.begin(); it != NULL; it = it->next, ++i) {
            a[it->key.first] = it->key;
            pos[it->key.first] = it->key.first;
        }
        cout << __LINE__ << endl;

        for (i = parent(size); i >= 0; --i) {
            minHeapify(i);
        }
    }

    ~MinHeap() {
        delete[] a;
        delete[] pos;
        a = NULL;
        size = 0;
    }

    int Size() {
        return size;
    }

    heapElem GetMinimum() {
        return a[0];
    }

    heapElem ExtractMin() {
        heapElem res = a[0];
        swap(pos[a[0].first], pos[a[size - 1].first]);
        a[0] = a[size - 1];
        --size;
        minHeapify(0);
        return res;
    }

    void DecreaseKey(const heapElem& key) {
        int i = pos[key.first];
        a[i] = key;
        while (i && cmpFun(a[i], a[parent(i)])) {
            swap(pos[a[i].first], pos[a[parent(i)].first]);
            swap(a[i], a[parent(i)]);
            i = parent(i);
        }
    }

    void InsertKey(const heapElem& key) {
        ++size;
        a[size - 1] = make_pair(key.first, 0);
        pos[key.first] = size - 1;
        DecreaseKey(key);
    }

    void DeleteKey(const heapElem& key) {
        int i = pos[key.first];
        printf("i = %d\n", i);
        swap(pos[a[i].first], pos[a[size - 1].first]);
        a[i] = a[size - 1];
        --size;
        minHeapify(i);
    }

private:
    void minHeapify(int i) {
        int l, r;
        l = left(i);
        r = right(i);

        int smallest;
        if (l < size && cmpFun(a[l], a[i])) {
            smallest = l;
        } else {
            smallest = i;
        }
        if (r < size && cmpFun(a[r], a[smallest])) {
            smallest = r;
        }

        if (smallest != i) {
            swap(pos[a[i].first], pos[a[smallest].first]);
            swap(a[i], a[smallest]);
            minHeapify(smallest);
        }
    }

    inline int left(int i) {
        return 1 + (i << 1);
    }

    inline int right(int i) {
        return 2 + (i << 1);
    }

    inline int parent(int i) {
        return ((i - 1) >> 1);
    }
};