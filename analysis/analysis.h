#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <core/common.h>
#include <core/vec.h>
#include <core/list.h>
#include <core/hash.h>
#include <core/paramvalue.h>

struct SolutionPoint {
	FP t;
	Vec<FP> y;

	bool operator<(const SolutionPoint& sp) const { return t < sp.t; }
};

struct Color {
	FP r, g, b;
};

struct SolutionLine {
	std::string _name;
	std::vector<SolutionPoint> _points;
	std::vector<Color> _colors;
	int _symbol;
};

struct ViewRect {
	double _xmin;
	double _xmax;
	double _ymin;
	double _ymax;
};

void ProcessAnalysisOptions(Hash<ParamValue>& params, List<ParamValue>& args);
void ReadFile(std::string filename, SolutionPoint& solutionPoint, Hash<ParamValue>& params);

extern std::vector<SolutionLine> g_SolutionLines;
extern std::vector<std::string> g_SolutionNames;
extern ViewRect g_ViewRect;

#endif

