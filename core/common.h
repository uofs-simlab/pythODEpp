#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <list>
#include <limits>

#include <cmath>
#include <cstring>
#include <cstdlib>

#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

typedef double FP;
typedef std::complex<FP> CFP;

inline void clamp(FP& v, FP mi, FP ma) {
	if( v < mi ) v = mi;
	else if( v > ma ) v = ma;
}

template <class T>
inline T sqr(T x) {
	return x*x;
}

inline FP fabs(CFP x) {
	return std::norm(x);
}

template <class T>
inline T realpow(T x, T y) {
	return pow(fabs(x), y) * cos( y*atan2((T)0,x) );
}

#endif

