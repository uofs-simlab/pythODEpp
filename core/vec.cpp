#include <core/vec.h>

Vec<FP> VecReal(const Vec<CFP>& cv) {
	Vec<FP> ret(cv.Size());
	for( long i = 0; i < cv.Size(); i++ )
		ret(i) = cv(i).real();
	return ret;
}

Vec<FP> VecImag(const Vec<CFP>& cv) {
	Vec<FP> ret(cv.Size());
	for( long i = 0; i < cv.Size(); i++ )
		ret(i) = cv(i).imag();
	return ret;
}

