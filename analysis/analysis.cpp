#include <core/exception.h>
#include <core/args.h>
#include <analysis/analysis.h>

std::string g_GraphTitle;
std::string g_GraphXLabel;
std::string g_GraphYLabel;
std::vector<SolutionLine> g_SolutionLines;
std::vector<std::string> g_SolutionNames;

Color g_ColorList[] = {
	{ 1.0, 0.0, 0.0 },
	{ 0.0, 0.7, 0.0 },
	{ 0.0, 0.0, 1.0 },
	{ 0.0, 0.9, 0.9 },
	{ 1.0, 0.0, 1.0 },
	{ 0.5, 0.3, 0.1 },
	{ 0.0, 0.0, 0.0 },
};
unsigned g_ColorCount = 7;

ViewRect g_ViewRect = { 0, 0, 0, 0 };
bool g_InitializedViewRect = false;

void ReadFile(std::string filename, SolutionPoint& solutionPoint, Hash<ParamValue>& params) {
	ParamValue* pv;
	long offset = 0; if( (pv = params.Get("solution offset")) ) offset = pv->GetLong();
	long stride = 1; if( (pv = params.Get("solution stride")) ) stride = pv->GetLong();
	long count  = -1; if( (pv = params.Get("solution count")) ) count = pv->GetLong();

	std::ifstream file;
	file.open(filename.c_str(), std::ios::binary);
	
	if( !file.is_open() )
		throw Exception() << "Unable to open " << filename << ".";
    
	file.read((char*)&solutionPoint.t, sizeof(solutionPoint.t));

	Vec<FP> tempVec;
	file >> tempVec;
	file.close();

	if( count < 0 ) count = ceil(FP(tempVec.Size()-offset)/stride);
	solutionPoint.y.Resize(count);

	for( long i = offset, j = 0; i < tempVec.Size() && j < count; i += stride, j++ ) {
		solutionPoint.y[j] = tempVec[i];
	}
}

std::string GetInputLine() {
	std::string line;
	if( !std::getline(std::cin, line) ) {
		if( std::cin.eof() )
			throw Exception() << "Input stream ended prematurely.";
		throw Exception() << "Unknown error in input stream.";
	}
	return line;
}

void LoadSolutions(Hash<ParamValue>& params) {
	if( !params.Get("number") )
		throw Exception() << "Number of solutions is required.";

	if( params.Get("solnames") ) {
		g_SolutionNames.resize(params["solnames"].GetLong());
		for( long i = 0; i < (long)g_SolutionNames.size(); i++ )
			g_SolutionNames[i] = GetInputLine();
	}

	long solsCount = params["number"].GetLong();
	g_SolutionLines.resize(solsCount);
	for( long s = 0, color = 0; s < solsCount; s++ ) {
		SolutionLine& sl = g_SolutionLines[s];

		std::string pathName = GetInputLine();
		sl._name = GetInputLine();
		sl._symbol = -1;

		std::vector<std::string> solPoints = ListDir(pathName.c_str());
		std::sort(solPoints.begin(), solPoints.end());

		sl._points.resize(solPoints.size());
		for( unsigned i = 0; i < solPoints.size(); i++ ) {
			SolutionPoint sp;
			ReadFile(solPoints[i], sp, params);
			sl._points[i] = sp;

			// Find margins
			if( !g_InitializedViewRect ) {
				g_ViewRect._xmin = sp.t;
				g_ViewRect._xmax = sp.t;
				g_ViewRect._ymin = sp.y[0];
				g_ViewRect._ymax = sp.y[0];
				g_InitializedViewRect = true;
			}

			if( sp.t < g_ViewRect._xmin )
				g_ViewRect._xmin = sp.t;
			if( sp.t > g_ViewRect._xmax )
				g_ViewRect._xmax = sp.t;
			for( long j = 0; j < sp.y.Size(); j++ ) {
				if( sp.y[j] < g_ViewRect._ymin )
					g_ViewRect._ymin = sp.y[j];
				if( sp.y[j] > g_ViewRect._ymax )
					g_ViewRect._ymax = sp.y[j];
			}
		}

		// Set a color for each point
		// These can be set the same for each point desired, takes just
		// a slight manipulation of the increment
		sl._colors.resize(sl._points[0].y.Size());
		for( unsigned c = 0; c < sl._colors.size(); c++, color++ ) {
			long ci = color%g_ColorCount;
			sl._colors[c].r = g_ColorList[ci].r;
			sl._colors[c].g = g_ColorList[ci].g;
			sl._colors[c].b = g_ColorList[ci].b;
		}
	}
}

void LoadSolution1D(Hash<ParamValue>& params) {
	std::string pathName = GetInputLine();
	std::vector<std::string> solPoints = ListDir(pathName.c_str());
	std::sort(solPoints.begin(), solPoints.end());

	long times = params["solution times"].GetLong();
	g_SolutionLines.resize(times);
	for( long t = 0; t < times; t++ ) {
		SolutionLine& sl = g_SolutionLines[t];

		sl._name = GetInputLine();
		sl._symbol = -1;
		
		sl._colors.resize(1);
		long ci = t%g_ColorCount;
		sl._colors[0].r = g_ColorList[ci].r;
		sl._colors[0].g = g_ColorList[ci].g;
		sl._colors[0].b = g_ColorList[ci].b;

		FP solTime;
		if( !(std::stringstream(sl._name) >> solTime) )
			throw Exception() << "Error reading time of solution.";
		sl._name = "time = " + sl._name;

		SolutionPoint spPrevious;
		ReadFile(solPoints[0], spPrevious, params);
		if( solTime < spPrevious.t )
			throw Exception() << "Time " << solTime << " is smaller than solution range beginning at " << spPrevious.t << ".";

		for( unsigned p = 1; ; p++ ) {
			if( p >= solPoints.size() )
				throw Exception() << "Requested time " << solTime << " is larger than the solution range ending at " << spPrevious.t << ".";

			SolutionPoint spCurrent;
			ReadFile(solPoints[p], spCurrent, params);

			// Fill solution line
			if( solTime <= spCurrent.t ) {
				FP alpha = (solTime-spPrevious.t)/(spCurrent.t-spPrevious.t);

				FP xmin = 0;
				FP xmax = (FP)(spCurrent.y.Size()-1);
				if( params.Get("xmin") ) xmin = params["xmin"].GetFP();
				if( params.Get("xmax") ) xmax = params["xmax"].GetFP();

				sl._points.resize(spCurrent.y.Size());
				for( unsigned xi = 0; xi < spCurrent.y.Size(); xi++ ) {
					FP x = xmin + xi*(xmax-xmin)/(spCurrent.y.Size()-1);
					FP y =  (1-alpha)*spPrevious.y[xi] + alpha*spCurrent.y[xi];
					sl._points[xi].t = x;
					sl._points[xi].y.Resize(1);
					sl._points[xi].y[0] = y;

					if( !g_InitializedViewRect ) {
						g_ViewRect._xmin = x;
						g_ViewRect._xmax = x;
						g_ViewRect._ymin = y;
						g_ViewRect._ymax = y;
						g_InitializedViewRect = true;
					}

					if( x < g_ViewRect._xmin )
						g_ViewRect._xmin = x;
					if( x > g_ViewRect._xmax )
						g_ViewRect._xmax = x;
					if( y < g_ViewRect._ymin )
						g_ViewRect._ymin = y;
					if( y > g_ViewRect._ymax )
						g_ViewRect._ymax = y;
				}
				break;
			}

			spPrevious = spCurrent;
		}
	}
}

void LoadSolution2D(Hash<ParamValue>& params) {
	std::string pathName = GetInputLine();
	std::vector<std::string> solPoints = ListDir(pathName.c_str());
	std::sort(solPoints.begin(), solPoints.end());

	g_SolutionLines.resize(1);
	FP solTime = params["solution time"].GetFP();
	SolutionLine& sl = g_SolutionLines[0];

	sl._symbol = -1;
	sl._points.resize(1);
		
	SolutionPoint spPrevious;
	ReadFile(solPoints[0], spPrevious, params);
	if( solTime < spPrevious.t )
		throw Exception() << "Time " << solTime << " is smaller than solution range beginning at " << spPrevious.t << ".";

	for( unsigned p = 1; ; p++ ) {
		if( p >= solPoints.size() )
			throw Exception() << "Requested time " << solTime << " is larger than the solution range ending at " << spPrevious.t << ".";

		SolutionPoint spCurrent;
		ReadFile(solPoints[p], spCurrent, params);

		if( solTime <= spCurrent.t ) {
			FP alpha = (solTime-spPrevious.t)/(spCurrent.t-spPrevious.t);
			sl._points[0].t = solTime;
			sl._points[0].y = (1-alpha)*spPrevious.y + alpha*spCurrent.y;
			break;
		}

		spPrevious = spCurrent;
	}
}

void LoadAccuracy(Hash<ParamValue>& params) {
	if( !params.Get("groups") )
		throw Exception() << "Number of groups is required.";
	long groupCount = params["groups"].GetLong();

	// Load reference solution
	SolutionPoint reference;	
	if( params.Get("reference") ) {
		std::vector<std::string> refFiles = ListDir(params["reference"].GetString());
		std::sort(refFiles.begin(), refFiles.end());
		ReadFile(refFiles.back(), reference, params);
	} else {
		long refSize;
		if( !(std::stringstream(GetInputLine()) >> refSize) )
			throw Exception() << "Error reading reference solution size.";
		// Read time
		if( !(std::stringstream(GetInputLine()) >> reference.t) )
			throw Exception() << "Error reading reference solution time.";

		reference.y.Resize(refSize);
		for( long i = 0; i < refSize; i++ ) {
			if( !(std::stringstream(GetInputLine()) >> reference.y[i]) )
				throw Exception() << "Error reading reference solution element.";
		}
	}

	g_SolutionLines.resize(groupCount);
	for( long g = 0, color = 0; g < groupCount; g++ ) {
		SolutionLine& sl = g_SolutionLines[g];
		long solNum;
		if( !(std::stringstream(GetInputLine()) >> solNum) )
			throw Exception() << "Error reading number of solutions.";
		sl._name = GetInputLine();
		sl._symbol = -1;
		if( params.Get("symbol") && params["symbol"].GetLong() )
			if( !(std::stringstream(GetInputLine()) >> sl._symbol) )
				throw Exception() << "Error reading solution symbol.";

		// Hack to allow for matching colors
		bool setColor = (bool)GetDefaultLong(params, "color", 0);
		FP setColorR, setColorG, setColorB;
		if( setColor )
			if( !(std::stringstream(GetInputLine()) >> setColorR >> setColorG >> setColorB) )
				throw Exception() << "Error reading solution color.";

		for( long s = 0; s < solNum; s++ ) {
			std::vector<std::string> solFiles = ListDir(GetInputLine().c_str());
			std::sort(solFiles.begin(), solFiles.end());

			// Store quantity
			SolutionPoint quantity;
			quantity.y.Resize(1);
			if( !(std::stringstream(GetInputLine()) >> quantity.y[0]) )
				throw Exception() << "Error loading quantity.";
			
			// Store accuracy
			SolutionPoint sp;
			ReadFile(solFiles.back(), sp, params);
			quantity.t = (sp.y-reference.y).Norm();

			// Discard NaNs
			if( quantity.t != quantity.t )
				continue;

			// Discard accuracy > 1
			if( quantity.t > 1 ) {
				//std::cout << "discarding accuracy (" << quantity.t << ") for " << sl._name << "\n";
				continue;
			}

			// Discard accuracy <= 0 because we can't plot it nicely on a log-log plot
			if( quantity.t <= 0 )
				continue;

			sl._points.push_back(quantity);

			// Find margins
			if( g_InitializedViewRect ) {
				if( quantity.t < g_ViewRect._xmin )
					g_ViewRect._xmin = quantity.t;
				if( quantity.t > g_ViewRect._xmax )
					g_ViewRect._xmax = quantity.t;
				if( quantity.y[0] < g_ViewRect._ymin )
					g_ViewRect._ymin = quantity.y[0];
				if( quantity.y[0] > g_ViewRect._ymax )
					g_ViewRect._ymax = quantity.y[0];
			} else {
					g_ViewRect._xmin = quantity.t;
					g_ViewRect._xmax = quantity.t;
					g_ViewRect._ymin = quantity.y[0];
					g_ViewRect._ymax = quantity.y[0];
					g_InitializedViewRect = true;
			}
		}
		
		// Set a color for each point
		if( sl._points.size() ) {
			sl._colors.resize(sl._points[0].y.Size());
			for( unsigned c = 0; c < sl._colors.size(); c++ ) {
				if( setColor ) {
					sl._colors[c].r = setColorR;
					sl._colors[c].g = setColorG;
					sl._colors[c].b = setColorB;
				} else {
					long ci = color%g_ColorCount;
					sl._colors[c].r = g_ColorList[ci].r;
					sl._colors[c].g = g_ColorList[ci].g;
					sl._colors[c].b = g_ColorList[ci].b;
					color++;
				}
			}
		}

		std::sort(sl._points.begin(), sl._points.end());
	}
}

void ProcessAnalysisOptions(Hash<ParamValue>& params, List<ParamValue>& args) {
	if( !params.Get("mode") )
		throw Exception() << "Analysis mode is required.";
	std::string modestr(params["mode"].GetString());

	if( modestr == "Solutions" )
		LoadSolutions(params);
	else if( modestr == "Solution1D" )
		LoadSolution1D(params);
	else if( modestr == "Solution2D" )
		LoadSolution2D(params);
	else if( modestr == "Portraits" )
		throw Exception() << "Phase portraits not yet implemented.";
	else if( modestr == "Accuracy" )
		LoadAccuracy(params);
	else
		throw Exception() << "Mode " << modestr << " is undefined.";
}


