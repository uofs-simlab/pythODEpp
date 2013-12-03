#ifndef TIMER_H
#define TIMER_H

#include <core/common.h>

class Timer {
	bool _running;
	timeval _startTime;
	timeval _endTime;

public:
	Timer() {
		Start();
	}

	void Start() {
		_running = true;
		gettimeofday(&_startTime, NULL);
	}

	void End() {
		_running = false;
		gettimeofday(&_endTime, NULL);
	}

	FP msec() {
		if( _running )
			gettimeofday(&_endTime, NULL);

		long seconds = _endTime.tv_sec-_startTime.tv_sec;
		long useconds = _endTime.tv_usec-_startTime.tv_usec;
		return seconds*1e3 + useconds/1e3;
	}

	FP sec() {
		if( _running )
			gettimeofday(&_endTime, NULL);

		long seconds = _endTime.tv_sec-_startTime.tv_sec;
		long useconds = _endTime.tv_usec-_startTime.tv_usec;
		return seconds + useconds/1e6;
	}
};

#endif

