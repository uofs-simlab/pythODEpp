#include <core/ioobject.h>

std::ostream& operator<<(std::ostream &out, const IOObject& obj) {
	obj.Dump(out);
	return out;
}

std::istream& operator>>(std::istream &in, IOObject& obj) {
	obj.Load(in);
	return in;
}

std::string IOObject::GenerateUID() {
	// Generate time part of the uid
	time_t rawtime;
	struct tm* timeinfo;
	char tbuf[80];
	time( &rawtime );
	timeinfo = localtime(&rawtime);
	strftime(tbuf, 80, "%Y-%m-%d--%H-%M-%S", timeinfo);
    
	// Calculate random part of the uid
	std::string uid;
    uid.resize(UID_LENGTH);
	
	static const char alphanum[] =
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    
	for( int i = 0; i < UID_LENGTH; i++ ) {
		uid[i] = alphanum[rand() % (sizeof(alphanum)-1)];
	}
    
	return std::string(tbuf) + "--" + uid;
}

void IOObject::MakePath(const char* path) {
    struct stat sb;
	
	char* end = (char*)strchr(path, '/');
	while( end ) {
		*end = '\0';
		if( stat(path, &sb) != 0 )
			mkdir(path, 504);
		*end = '/';
		end = strchr(end+1, '/');
	}
	
    if( stat(path, &sb) != 0 )
        mkdir(path, 504);
}

