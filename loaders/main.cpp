#include <core/common.h>
#include <core/exception.h>
#include <core/args.h>
#include <analysis/analysis.h>
#include <ivps/baseivp.h>
#include <methods/basemethod.h>
#include <solvers/basesolver.h>

BaseIVP* AllocIVP(Hash<ParamValue>& params);
BaseMethod* AllocMethod(Hash<ParamValue>& params, BaseIVP* ivp);
BaseSolver* AllocSolver(Hash<ParamValue>& params, BaseMethod* method, BaseIVP* ivp);

void EigenvalueMain(Hash<ParamValue>& params, List<ParamValue>& args) {
	// Get the solution points
	if( !params.Get("path") )
		throw Exception() << "Path is required.";
	std::vector<std::string> solPoints = ListDir(params["path"].GetString());
	std::sort(solPoints.begin(), solPoints.end());
	std::reverse(solPoints.begin(), solPoints.end());

	// Get the IVP
	if( !params.Get("ivp") )
		throw Exception() << "IVP class is required.";

	// Allocate the IVP
	BaseIVP* ivp;
	if( !(ivp = AllocIVP(params)) )
		throw Exception() << "IVP class " << params["ivp"].GetString() << " is not defined.";
	ivp->InitializeDerivatives();

	if( !params.Get("dt") )
		throw Exception() << "Timestep dt is required.";
	FP dt = params["dt"].GetFP();

	SolutionPoint spPrev;
	ReadFile(solPoints.back(), spPrev, params);
	solPoints.pop_back();

	// Output python header
	std::cout << "from scipy.sparse import *\n";
	std::cout << "from scipy.linalg import *\n";
	std::cout << "from scipy import *\n";
	std::cout << "eigs = []\n";

	FP t = 0;
	while( solPoints.size() ) {
		SolutionPoint spCurrent;
		ReadFile(solPoints.back(), spCurrent, params);
		solPoints.pop_back();

		if( spCurrent.t > t ) {
			FP alpha = (t-spPrev.t)/(spCurrent.t-spPrev.t);
			CSRMat<FP>* jac = (CSRMat<FP>*)ivp->JacSparse(t, (1-alpha)*spPrev.y + alpha*spCurrent.y);
			jac->PrintPython(std::cout);

			std::cout << "eigs = eigs + [e for e in eig(csr_matrix((array(s),array(j),array(i)), shape=(m,n)).todense())[0] ]\n";
			t += dt;
		}

		spPrev = spCurrent;
	}

	std::cout << "print \"unset key\"\n";
	std::cout << "print \"set terminal pdf\"\n";
	std::cout << "print \"plot '-'\"\n";
	std::cout << "for e in eigs:\n";
	std::cout << "  print e.real, e.imag\n";
	std::cout << "print 'e'\n";

	delete ivp;
}

// -----------------------------------------------------------------------------

void DumpSolutionMain(Hash<ParamValue>& params, List<ParamValue>& args) {
	if( !params.Get("path") )
		throw Exception() << "Path is required.";

	std::vector<std::string> solPoints = ListDir(params["path"].GetString());
	std::sort(solPoints.begin(), solPoints.end());

	SolutionPoint sp;
	ReadFile(solPoints.back(), sp, params);
	std::cout << "t = " << sp.t << "\n";
	sp.y.Print(std::cout);
}

// -----------------------------------------------------------------------------

void RunnerMain(Hash<ParamValue>& params, List<ParamValue>& args) {
	BaseIVP* ivp = 0;
	BaseMethod* method = 0;
	BaseSolver* solver = 0;

	if( !params.Get("ivp") )
		throw Exception() << "IVP class is required.";

	if( !params.Get("method") )
		throw Exception() << "method class is required.";

	if( !params.Get("solver") )
		throw Exception() << "solver class is required.";

	if( !(ivp = AllocIVP(params)) )
		throw Exception() << "IVP class " << params["ivp"].GetString() << " is not defined.";
	ivp->InitializeDerivatives();

	if( !(method = AllocMethod(params, ivp)) )
		throw Exception() << "method class " << params["method"].GetString() << " is not defined.";

	if( !(solver = AllocSolver(params, method, ivp)) )
		throw Exception() << "solver class " << params["solver"].GetString() << " is not defined.";

	solver->RunSimulation();
	solver->DumpRunInfo(params);
}

// -----------------------------------------------------------------------------

#include <methods/rk.h>

void DumpTableauMain(Hash<ParamValue>& params, List<ParamValue>& args) {
	bool ark = false;

	// Better hope it is an RK method
	RKMethod* method = (RKMethod*)AllocMethod(params, 0);

	const char* methodName = params["ivp"].GetString();
	if( strcmp(methodName, "ARK3") == 0 ||
		strcmp(methodName, "ARK4") == 0 ||
		strcmp(methodName, "ARK5") == 0 )
		ark = true;

	for( ListNode<ParamValue>* pvn = args.Head(); pvn; pvn = pvn->_next ) {
		if( std::string((**pvn).GetString()) == "A" ) {
			method->GetA().Print(std::cout);
		} else if( std::string((**pvn).GetString()) == "b" ) {
			method->GetB().Print(std::cout);
		} else if( std::string((**pvn).GetString()) == "c" ) {
			method->GetC().Print(std::cout);
		} else if( std::string((**pvn).GetString()) == "A2" ) {
			if( ark )
				((IMEX*)method)->GetA2().Print(std::cout);
			else
				method->GetA().Print(std::cout);
		} else if( std::string((**pvn).GetString()) == "b2" ) {
			if( ark )
				((IMEX*)method)->GetB2().Print(std::cout);
			else
				method->GetB().Print(std::cout);
		} else if( std::string((**pvn).GetString()) == "c2" ) {
			if( ark )
				((IMEX*)method)->GetC2().Print(std::cout);
			else
				method->GetC().Print(std::cout);
		} else {
			throw Exception() << "Unrecognized coefficient '" << (**pvn).GetString() << "'.\n";
		}
	}

	delete method;
}

// -----------------------------------------------------------------------------

void GNUPlotMain(Hash<ParamValue>& params, List<ParamValue>& args);

int main(int argc, char* argv[]) {
	srand(time(0));
	int returnCode = 0;
	
	try {
		Hash<ParamValue> params;
		List<ParamValue> args;
		ProcessArgs(params, args, argc, argv);
	
		std::cout << std::fixed << std::setprecision(GetDefaultLong(params,"precision",30));

		if( !params.Get("phase") )
			throw Exception() << "Simulation phase is required.";

		const char* phase = params["phase"].GetString();

		if( strcmp(phase, "runner") == 0 )
			RunnerMain(params, args);
		else if( strcmp(phase, "gnuplot1d") == 0 )
			GNUPlotMain(params, args);
		else if( strcmp(phase, "eigenvalues") == 0 )
			EigenvalueMain(params, args);
		else if( strcmp(phase, "dumptableau") == 0 )
			DumpTableauMain(params, args);
		else if( strcmp(phase, "dumpsolution") == 0 )
			DumpSolutionMain(params, args);
		else
			throw Exception() << "Unrecognized phase " << phase << ".";
	} catch(Exception e) {
		std::cerr << (const char*)e << std::endl;
		returnCode = 1;
	}

#ifdef DEBUGBUILD
	Vec<FP>::PrintStats();
	Mat<FP>::PrintStats();
#endif

	return returnCode;
}
