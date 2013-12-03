#ifndef HASH_H
#define HASH_H

template <class T>
class Hash {
	class Bin {
	public:
		Bin(const char* name) {
			int l = 1 + strlen(name);
			_name = new char[l];
			strcpy(_name, name);
			_next = 0;
		}
		~Bin() { delete [] _name; }
		char* _name;
		T _obj;
		
		Bin* _next;
	};
	
	Bin** _bins;
	long _binCount;
	void (*_deleteCallback)(T&);
	
	unsigned GetBin(const char* key) const {
		unsigned long index = 0;
		while( *key ) {
			index = 33*index + *key;
			key++;
		}
		
		return index % _binCount;
	}
	
public:
	class Iterator {
		long _currentBin;
		Bin* _bin;
		const Hash* _h;
		
	public:
		Iterator(const Hash& h) : _h(&h) {
			_currentBin = 0;
			_bin = 0;
			while( _currentBin < _h->_binCount ) {
				if( (_bin = _h->_bins[_currentBin]) )
					break;
				_currentBin++;
			}
		}
		
		void Next() {
			if( _bin && _bin->_next) {
				_bin = _bin->_next;
				return;
			}
			
			_bin = 0;
			while( _currentBin < _h->_binCount-1 ) {
				_currentBin++;
				if( (_bin = _h->_bins[_currentBin]) )
					break;
			}
		}
		
		T* operator*() {
			if( _bin )
				return &_bin->_obj;
			return 0;
		}
		
		const char* CurrentKey() {
			if( _bin )
				return _bin->_name;
			return 0;
		}
	};
	
	Hash(void (*callback)(T&)=0, long binCount=10) {
		_binCount = binCount;
		_bins = new Bin*[_binCount];
		_deleteCallback = callback;
		memset(_bins, 0, _binCount*sizeof(Bin*));
	}
	
	Hash(const Hash& h) : _binCount(h._binCount) {
		_bins = new Bin*[_binCount];
		_deleteCallback = h._deleteCallback;
		memset(_bins, 0, _binCount*sizeof(Bin*));
		
		Iterator it(h);
		while( *it ) {
			(*this)[it.CurrentKey()] = **it;
			it.Next();
		}
	}
	
	~Hash() {
		Empty();
		if( _bins )
			delete [] _bins;
	}
	
	void Empty() {
		for( long i = 0; i < _binCount; i++ ) {
			Bin* b = _bins[i];
			
			while( b ) {
				Bin* t = b;
				if( _deleteCallback )
					_deleteCallback(b->_obj);
				b = b->_next;
				delete t;
			}
			
			_bins[i] = 0;
		}
	}
	
	T* Get(const char* key) {
		if( !key ) return 0;
		Bin* b = _bins[GetBin(key)];

		while( b ) {
			if( strcmp(b->_name,key) == 0 )
				return &(b->_obj);
			b = b->_next;
		}
		
		return 0;
	}
	
	T& operator [](const char* key) {
		int i = GetBin(key);
		Bin* b = _bins[i];
		while( b ) {
			if( strcmp(b->_name,key) == 0 )
				return b->_obj;
			b = b->_next;
		}
	
		b = new Bin(key);
		b->_next = _bins[i];
		_bins[i] = b;
		return b->_obj;
	}
	
	void Remove(const char* key) {
		int i = GetBin(key);
		
		Bin** b = _bins + i;
		
		while( *b ) {
			if( strcmp((*b)->_name, key) == 0 ) {
				Bin* t = *b;
				(*b) = (*b)->_next;
				
				if( _deleteCallback )
					_deleteCallback(t->_obj);
				delete t;
				return;
			}
			
			b = &(*b)->_next;
		}
	}
	
	void SetDeleteCallback(void (*callback)(T&)) { _deleteCallback = callback; }
};

#endif
