#ifndef LOADER_H
#define LOADER_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <boost/filesystem/operations.hpp>
#include "computation/module.H"
#include "computation/optimization/simplifier.H" // for simplifier_options

class module_loader: public simplifier_options
{
    mutable std::map<std::string, Module> modules;

public:
    std::vector<std::string> plugins_path;

    bool try_add_plugin_path(const std::string& path);
  
    std::string find_module(const std::string& modid) const;
    std::string find_plugin(const std::string& plugin_name) const;
    /// Load the module file from disk.
    /// Add the hard-coded elements to the module.
    Module load_module_from_file(const std::string&) const;
    Module load_module(const std::string&) const;
    module_loader(const std::vector<boost::filesystem::path>& paths);
};

expression_ref load_builtin(const module_loader& L, const std::string& symbol_name, const std::string& filename, int n, const std::string& fname);
expression_ref load_builtin(const std::string& symbol_name, const std::string& filename, int n, const std::string& fname);

extern const std::string plugin_extension;
#endif
