#include <core/exception.h>
#include <analysis/analysis.h>
#include <core/args.h>

#include <stdarg.h>

FILE* g_pdfFP = 0;
FILE* g_textFP = 0;

std::string RGBToHex(FP r, FP g, FP b) {
	clamp(r,0,1);
	clamp(g,0,1);
	clamp(b,0,1);
	long rl = long(255*r + 0.5);
	long gl = long(255*g + 0.5);
	long bl = long(255*b + 0.5);
	long rgb = ((rl & 0xff) << 16) + ((gl & 0xff) << 8) + (bl & 0xff);

	std::stringstream ss;
	ss << std::hex << std::setfill('0') << std::setw(6) << rgb;
	return ss.str();
}

void EchoPrint(const char* msg, ...) {
	va_list ptr;

	if( g_pdfFP ) {
		va_start(ptr, msg);
		vfprintf(g_pdfFP, msg, ptr);
		va_end(ptr);
	}

	if( g_textFP ) {
		va_start(ptr, msg);
		vfprintf(g_textFP, msg, ptr);
		va_end(ptr);
	}
	
	if( 0 ) {
		va_start(ptr, msg);
		vprintf(msg, ptr);
		va_end(ptr);
	}
}

void SetupPlot(Hash<ParamValue>& params) {
	FP xsize = 11, ysize = 8.5;
	if( params.Get("xsize") ) {
		xsize = params["xsize"].GetFP();
		ysize = xsize*8.5/11;
	}
	if( params.Get("ysize") ) {
		ysize = params["ysize"].GetFP();
		if( !params.Get("xsize") )
			xsize = ysize*11/8.5;
	}

	EchoPrint("set terminal pdf size %.2f,%.2f\n", xsize, ysize);
	EchoPrint("set termopt enhanced\n");
	EchoPrint("set output '%s'\n", params["plotfile"].GetString());
	EchoPrint("set key t rm\n");
}

void CUSPPlot(Hash<ParamValue>& params) {
	if( g_SolutionLines.size() != 1 )
		throw Exception() << "CUSP plots only makes sense with 1 solution vector.\n";

	SetupPlot(params);
	
	SolutionLine& sl = g_SolutionLines[0];
	if( sl._points[0].y.Size() % 3 )
		throw Exception() << "CUSP plots must have solution size divisible by 3. Size is " << sl._points[0].y.Size() << "\n";

	EchoPrint("unset key\n");
	EchoPrint("set lmargin at screen 0.1\n");
	EchoPrint("set rmargin at screen 0.9\n");
	EchoPrint("set tmargin at screen 0.9\n");
	EchoPrint("set bmargin at screen 0.1\n");
	EchoPrint("set yzeroaxis lt -1 lw 3\n");
	EchoPrint("set xzeroaxis lt -1 lw 3\n");
	EchoPrint("set zzeroaxis lt -1 lw 3 \n");
	EchoPrint("set xyplane at 0\n");
	EchoPrint("set format x \"\"\n");
	EchoPrint("set format y \"\"\n");
	EchoPrint("set xtics 0.25\n");
	EchoPrint("set ytics 0.25\n");
	EchoPrint("set xrange [-3:3]\n");
	EchoPrint("set yrange [-3:3]\n");
	EchoPrint("set zrange [-3:3]\n");
	EchoPrint("set grid\n");
	EchoPrint("set view 60,100\n");
	EchoPrint("unset ztics\n");
	EchoPrint("unset border\n");

	EchoPrint("splot ");
	bool comma = false;
	for( long s = 0; s < sl._points[0].y.Size(); s++ ) {
		if( comma ) EchoPrint(", ");
		else comma = true;

		EchoPrint("\"-\" with lines lw 3 lc rgb '#000000'");
	}
	EchoPrint("\n");

	for( long s = 0; s < sl._points[0].y.Size(); s += 3 ) {
		for( long l = 0; l < sl._points.size(); l++ )
			EchoPrint("%g %g %g\n", sl._points[l].y[s+1], sl._points[l].y[s+2], sl._points[l].y[s]);
		EchoPrint("e\n");
	}
}

void Plot2D(Hash<ParamValue>& params) {
	EchoPrint("set grid\n");
	EchoPrint("set border linewidth 1.5\n");
	if( params["logscale"].GetLong() )
		EchoPrint("set logscale xy\n");
	if( !params.Get("linetype") )
		params["linetype"].SetString("linespoints");

	// Write out linestyles
	unsigned lineCounter = 0;
	for( unsigned l = 0, ls = 1; l < g_SolutionLines.size(); l++ ) {
		SolutionLine& sl = g_SolutionLines[l];
		for( unsigned c = 0; c < sl._colors.size(); c++, ls++ ) {
			lineCounter++;
			std::string colorstr = RGBToHex(sl._colors[c].r, sl._colors[c].g, sl._colors[c].b);
			EchoPrint("set style line %d lc rgb '#%s' lt 1 lw 2 pt %d ps 1.5\n", ls, colorstr.c_str(), sl._symbol < 0 ? ls : sl._symbol);
		}
	}

	if( lineCounter == 0 )
		throw Exception() << "No solutions will be plotted.";
	SetupPlot(params);
	
	// Initiate plot command
	EchoPrint("plot ");
	bool comma = false;
	for( unsigned l = 0, ls = 1; l < g_SolutionLines.size(); l++ ) {
		SolutionLine& sl = g_SolutionLines[l];

		if( !sl._points.size() )
			continue;

		for( long e = 0; e < sl._points[0].y.Size(); e++, ls++ ) {
			if( comma ) EchoPrint(", ");
			else comma = true;

			std::string title = sl._name;
			if( g_SolutionNames.size() )
				title = g_SolutionNames[e] + " " + title;

			EchoPrint("\"-\" with %s ls %d title \"%s\"", params["linetype"].GetString(), ls, title.c_str());
		}
	}
	EchoPrint("\n");

	// Write them all out
	for( unsigned l = 0; l < g_SolutionLines.size(); l++ ) {
		SolutionLine& sl = g_SolutionLines[l];

		if( !sl._points.size() )
			continue;

		for( long e = 0; e < sl._points[0].y.Size(); e++ ) {
			for( unsigned j = 0; j < sl._points.size(); j++ )
				EchoPrint("%g %g\n", sl._points[j].t, sl._points[j].y[e]);
			EchoPrint("e\n");
		}
	}
}

void PlotSurface(Hash<ParamValue>& params) {
	if( g_SolutionLines.size() != 1 )
		throw Exception() << "One and only one solution is supported.";
	SetupPlot(params);

	//EchoPrint("set dgrid3d 50,50\n");
	EchoPrint("set hidden3d\n");
	EchoPrint("unset key\n");

	SolutionLine& sl = g_SolutionLines[0];
	EchoPrint("set xrange [%.3f:%.3f] reverse\n", sl._points[0].t, sl._points[sl._points.size()-1].t);

	// Setup units for the solution vector
	FP ymin = 0;
	FP ymax = sl._points[0].y.Size();
	if( params.Get("ymin") ) ymin = params["ymin"].GetFP();
	if( params.Get("ymax") ) ymax = params["ymax"].GetFP();
	FP dy = (ymax-ymin)/sl._points[0].y.Size();

	if( params.Get("ytics") )
		EchoPrint("set ytics %.4f\n", params["ytics"].GetFP());

	EchoPrint("splot '-' with lines\n");
	for( long t = 0; t < (long)sl._points.size(); t++ ) {
		SolutionPoint& sp = sl._points[t];
		for( long x = 0; x < sp.y.Size(); x++ ) {
			EchoPrint("%.3f %.3f %.3f\n", sp.t, dy*x, sp.y[x]);
		}
		EchoPrint("\n");
	}
	EchoPrint("e\n");
}

void ColorPlot(Hash<ParamValue>& params) {
	SetupPlot(params);

	SolutionLine& sl = g_SolutionLines[0];
	SolutionPoint& sp = sl._points[0];

	if( params.Get("cbmin") && params.Get("cbmax") )
		EchoPrint("set cbrange [%.3f:%.3f]\n", params["cbmin"].GetFP(), params["cbmax"].GetFP());

	long xdim = params["xdim"].GetLong();
	long ydim = params["ydim"].GetLong();

	EchoPrint("set xrange [%.3f:%.3f]\n", 0., (FP)(xdim-1));
	EchoPrint("set yrange [%.3f:%.3f]\n", 0., (FP)(ydim-1));
	EchoPrint("unset key\n");

	EchoPrint("plot '-' with image\n");
	for( long x = 0; x < xdim; x++ ) {
		for( long y = 0; y < ydim; y++ ) {
			EchoPrint("%.3f %.3f %.3f\n", (FP)x, (FP)y, sp.y[y*xdim+x]);	
		}
	}
	EchoPrint("e\n");
}

void GNUPlotMain(Hash<ParamValue>& params, List<ParamValue>& args) {
	ProcessAnalysisOptions(params, args);

	if( !params.Get("plotfile") )
		throw Exception() << "plotfile is required.";

	std::cout << "Plotting " << params["plotfile"].GetString() << "." << std::endl;

	// Open handle to gnuplot	
	if( GetDefaultLong(params, "plotpdf", 1) ) {
		if( !(g_pdfFP = popen("gnuplot", "w")) )
			throw Exception() << "Unable to open gnuplot.";
	}
	if( GetDefaultLong(params, "plottxt", 0) ) {
		if( !(g_textFP = fopen((std::string(params["plotfile"].GetString()) + ".txt").c_str(), "w")) )
			throw Exception() << "Unable to open plot text file.";
	}

	EchoPrint("reset\n");

	ParamValue* pv;
	if( (pv = params.Get("title")) ) EchoPrint("set title \"%s\"\n", pv->GetString());
	if( (pv = params.Get("xlabel")) ) EchoPrint("set xlabel \"%s\"\n", pv->GetString());
	if( (pv = params.Get("ylabel")) ) EchoPrint("set ylabel \"%s\"\n", pv->GetString());
	if( (pv = params.Get("zlabel")) ) EchoPrint("set zlabel \"%s\"\n", pv->GetString());

	if( params["surface"].GetLong() )
		PlotSurface(params);
	else if( params["colormap"].GetLong() )
		ColorPlot(params);
	else if( params["cuspplot"].GetLong() )
		CUSPPlot(params);
	else
		Plot2D(params);

	if( g_pdfFP ) pclose(g_pdfFP);
	if( g_textFP ) fclose(g_textFP);
}

