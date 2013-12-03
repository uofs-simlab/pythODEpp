#ifndef IO_OBJECT_H
#define IO_OBJECT_H

#include <core/common.h>

#define UID_LENGTH 8

class IOObject {
protected:
	virtual void Dump(std::ostream &out) const = 0;
	virtual void Load(std::istream &in) = 0;
	
public:
	static std::string GenerateUID();
	static void MakePath(const char* path);
	
	friend std::ostream& operator<<(std::ostream& out, const IOObject& obj);
	friend std::istream& operator>>(std::istream& in, IOObject& obj);
};

std::ostream& operator<<(std::ostream& out, const IOObject& obj);
std::istream& operator>>(std::istream& in, IOObject& obj);

#endif
