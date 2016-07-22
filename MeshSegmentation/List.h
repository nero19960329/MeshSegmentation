#pragma once

template <class listElem>
class ListNode {
public:
    listElem key;
    ListNode *next;

public:
    ListNode() : next(NULL) {}
    ListNode(const listElem& elem) : key(elem), next(NULL) {}
};

template <class listElem>
class List {
private:
    ListNode<listElem> *root, *tail;
    int size;

public:
    typedef ListNode<listElem>* iterator;

public:
    List() : size(0) {
        root = new ListNode<listElem>();
        tail = root;
    }

    void push_back(const listElem& elem) {
        ListNode<listElem> *tmp = new ListNode<listElem>(elem);
        tail->next = tmp;
        tail = tmp;
        ++size;
    }

    int Size() { return size; }
    iterator begin() { return root->next; }
};