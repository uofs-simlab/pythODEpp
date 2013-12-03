#ifndef LIST_H
#define LIST_H

#include <stdint.h>

template <class T>
class ListNode {
public:
	ListNode(const T& obj) : _next(0), _prev(0), _object(obj) { }
	ListNode<T> *_next, *_prev;
	T _object;
    
    T& operator*() { return _object; }
};

template <class T>
class List {
	ListNode<T> *_head, *_tail;
	long _count;
	void (*_deleteCallback)(T&);
	
	List(const List& list) { }
    
public:
	List() : _head(0), _tail(0), _count(0), _deleteCallback(0) { }
	List(void (*deleteCallback)(T&)) : _head(0), _tail(0), _count(0),_deleteCallback(deleteCallback) { }
	~List() { Empty(); }
    
	T Pop() {
		T obj = _tail->_object;
		
		if( _head == _tail ) {
			delete _head;
			_head = _tail = 0;
			_count = 0;
		} else {
			_tail = _tail->_prev;
			delete _tail->_next;
			_tail->_next = 0;
			_count--;
		}
		
		return obj;
	}
	
	T Shift() {
		T obj = _head->_object;
		
		if( _head == _tail ) {
			delete _head;
			_head = _tail = 0;
			_count = 0;
		} else {
			_head = _head->_next;
			delete _head->_prev;
			_head->_prev = 0;
			_count--;
		}
		
		return obj;
	}
    
	void Push(const T& obj) {
		if(!_head) {
			_head = new ListNode<T>(obj);
			_tail = _head;
			_count = 1;
		} else {
			_tail->_next = new ListNode<T>(obj);
			_tail->_next->_prev = _tail;
			_tail = _tail->_next;
			_count++;
		}
	}
	
	void Unshift(const T& obj) {
		if(!_head) {
			_head = new ListNode<T>(obj);
			_tail = _head;
			_count = 1;
		} else {
			_head->_prev = new ListNode<T>(obj);
			_head->_prev->_next = _head;
			_head = _head->_prev;
			_count++;
		}
	}
	
	void InsertBefore(ListNode<T>* curNode, const T& newObj) {
		if( !curNode ) {
			Push(newObj);
			return;
		}
        
		if( _head == curNode ) {
			Unshift(newObj);
			return;
		}
        
		ListNode<T>* newNode = new ListNode<T>(newObj);
		newNode->_next = curNode;
		newNode->_prev = curNode->_prev;
		curNode->_prev = newNode;
		newNode->_prev->_next = newNode;
	}
    
	T Remove(ListNode<T>* node) {
		if( _head == node )
			return Shift();
		if( _tail == node )
			return Pop();
        
		T obj = node->_object;
		node->_next->_prev = node->_prev;
		node->_prev->_next = node->_next;
		delete node;
        
		return obj;
	}
    
	void Empty() {
		while( _head ) {
			if( _deleteCallback )
				_deleteCallback(_head->_object);
			ListNode<T>* temp = _head;
			_head = _head->_next;
			delete temp;
		}
		_count = 0;
		_tail = 0;
	}
    

public:
	void SetDeleteCallback(void (*callback)(T&)) { _deleteCallback = callback; }
    
	ListNode<T>* Head() { return _head; }
	ListNode<T>* Tail() { return _tail; }
	long Count() { return _count; }
};

#endif

