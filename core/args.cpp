#include <core/exception.h>
#include <core/args.h>

#include <dirent.h>

void ProcessArgs(Hash<ParamValue>& params, List<ParamValue>& args, int argc, char* argv[]) {
	for( int a = 1; a < argc; a++ ) {
		if( !argv[a][0] )
			continue;
		
		if( argv[a][0] == '-' ) { 
			const char* key = argv[a]+1;

			if( !*key )
				continue;
			
			if( a+1 >= argc )
				throw Exception() << key << " requires a parameter.";
			
			params[key].SetString(argv[++a]);
		} else {
			args.Push(ParamValue(argv[a]));
		}
	}
}

std::vector<std::string> ListDir(const char* path) {
	std::vector<std::string> filelist;
	
	DIR* d = opendir(path);
	if( !d )
		throw Exception() << "Unable to open path " << path << ".";

	struct dirent* dir;
	while( (dir = readdir(d)) )
		if( dir->d_name[0] != '.' )
			filelist.push_back(path + ("/" + std::string(dir->d_name)));
	closedir(d);
	
	return filelist;
}

