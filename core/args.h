#ifndef ARGS_H
#define ARGS_H

#include <core/hash.h>
#include <core/list.h>
#include <core/paramvalue.h>

void ProcessArgs(Hash<ParamValue>& params, List<ParamValue>& args, int argc, char* argv[]);
std::vector<std::string> ListDir(const char* path);

#endif

