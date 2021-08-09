

#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <filesystem>
#include "ParamBin/parambin.hpp"
#include "ParamBin/parambinstrings.h"

template<>
void ParamBin::get(const std::string& name,double& val) const
{
    std::string strval = getStrParam(name);
    std::string rawval = pwbin::eatWhiteSpace(pwbin::removeSubString(strval,'[',']'));
    convertFromString<double>(rawval,val);
}

template<>
void ParamBin::get(const std::string& name,std::vector<double>& vals) const
{
    std::string strval = getStrParam(name);
    std::vector<std::string> strvec = pwbin::parseString(strval,',');
    vals.clear();
    for(auto strval : strvec){
        double val;
        std::string rawval = pwbin::eatWhiteSpace(pwbin::removeSubString(strval,'[',']'));
        convertFromString<double>(rawval,val);
        vals.push_back(val);
    }
}

// Only children are copied, parent bin is considered root of the tree
ParamBin::ParamBin(const ParamBin& bin) : m_parent(nullptr),
    m_params(bin.m_params),
    m_children(),
    m_alias_map(bin.m_alias_map),
    m_ralias_map(bin.m_ralias_map)
{
    for(auto it = bin.m_children.cbegin() ; it != bin.m_children.cend(); it++)
    {
        std::string bin_name = it->first;
        std::unique_ptr<ParamBin> child_bin(new ParamBin(*it->second));
        child_bin->m_parent = this;
        m_children[bin_name] = std::move(child_bin);
    }
}

// Assignment
ParamBin& ParamBin::operator = (const ParamBin& bin)
{
    m_parent = nullptr;
    m_params = bin.m_params;
    m_alias_map =  bin.m_alias_map;
    m_ralias_map = bin.m_ralias_map;
    {
        for(auto it = bin.m_children.cbegin() ; it != bin.m_children.cend(); it++)
        {
            std::string bin_name = it->first;
            std::unique_ptr<ParamBin> child_bin(new ParamBin(*it->second));
            child_bin->m_parent = this;
            m_children[bin_name] = std::move(child_bin);
        }
    }
    return *this;
}

ParamBin::ParamBin() : 
    m_parent(nullptr) {}

ParamBin::ParamBin(const char* FILE) : 
    m_parent(nullptr) { loadParamFile(FILE); }

ParamBin::ParamBin(std::string fileName) : 
    m_parent(nullptr) { loadParamFile(fileName); }

ParamBin& ParamBin::operator<<(const NamedBin& named_bin)
{
    setBin(named_bin.name(),named_bin.bin());
    return *this;
}

std::string ParamBin::getStrParam(const std::string& name) const
{
    std::string strval;
    std::string key = pwbin::eatWhiteSpace(name);
    if(searchParamMap(key,strval))
        return strval;
    if(rootAliasSearch(key,strval)) 
        return strval;
    throw ParamBinKeyException(key);
}

void ParamBin::set(const NamedBin& named_bin) 
{
    setBin(named_bin.name(),named_bin.bin());
}

void ParamBin::setAlias(const std::string& name,const std::string& alias)
{
    std::string key = pwbin::eatWhiteSpace(name);
    if(m_params.count(key) > 0){
        m_alias_map[alias] = key;
        m_ralias_map[key] = alias;
    } else
        throw ParamBinKeyException(key);
}

std::string ParamBin::setParamKey(const std::string& name)
{
    std::string key = pwbin::eatWhiteSpace(pwbin::removeSubString(name,'(',')'));
    std::string alias = pwbin::eatWhiteSpace(pwbin::subString(name,'(',')')); 
    if(!alias.empty()){
        m_alias_map[alias] = key;
        m_ralias_map[key] = alias;
    }
    return key;
}

bool ParamBin::searchParamMap(const std::string& key,std::string& strval) const
{
    auto it = m_params.find(key);
    if(it != m_params.cend()){
        strval = (*it).second;
        return true;
    }
    return false;
}

bool ParamBin::searchAliasTree(const std::string& alias_key,std::string& strval) const
{
    // Check current bin alias'
    auto it = m_alias_map.find(alias_key);
    if(it != m_alias_map.cend()){
        std::string key = (*it).second;
        return searchParamMap(key,strval);
    }
    // Check child bins
    for(const auto& kv : m_children) {
        if(kv.second->searchAliasTree(alias_key,strval))
            return true;
    }
    return false;
}

bool ParamBin::rootAliasSearch(const std::string& alias_key,std::string& strval) const 
{
    const ParamBin* top = this;
    while(top->m_parent != nullptr)
        top = top->m_parent;
    return top->searchAliasTree(alias_key,strval);
}


template <>
void convertFromString(const std::string& str,char& val) { val = *str.c_str(); }

template <>
void convertFromString(const std::string& str,std::string& val) { val = str; }

template <>
void convertFromString(const std::string& str,std::vector<std::string>& vals)
{
    vals.clear();
    vals = pwbin::parseString(str,',');
}

template <>
std::string convertToString(std::string str) { return str; }

template <>
std::string convertToString(const char* cstr) { return std::string(cstr); }

template <>
std::string convertToString(char val) { return std::string(1,val); }



void lineToNameVal(const std::string& line_feed,std::string& name,std::string& vals)
{
    std::vector<std::string> parsedParam = pwbin::parseString(line_feed,':');
    name = parsedParam[0];
    vals = parsedParam[1];
    name = pwbin::eatWhiteSpace(name);
    vals = pwbin::eatWhiteSpace(vals);
}

std::ifstream& readNextLine(std::ifstream& fin,std::string& line_feed) 
{
    getline(fin,line_feed);
    line_feed = pwbin::decommentString(line_feed,"#");
    // Check that line_feed is not whitespace
    while(fin && pwbin::isWhitespace(line_feed)){
        getline(fin,line_feed);
        line_feed = pwbin::decommentString(line_feed,"#");
    }

    return fin;
}

void ParamBin::scanYAML(std::ifstream& fin)
{
    // Levels determine parent/child relationships, root level is -1 with a parent=nullptr
    int level = -1;
    std::vector<int> levels;
    std::vector<ParamBin*> parents;
    levels.push_back(level);

    parents.push_back(this);
		std::string line_feed;

    while(readNextLine(fin,line_feed))
    {
				std::vector<std::string> parsed_line_feed = pwbin::parseString(line_feed,':');
				if(parsed_line_feed.size() == 0)
						throw std::exception(); 
				else if(parsed_line_feed.size() == 1){
						throw std::exception(); 
				} else if(parsed_line_feed.size() == 2){
						std::string name(pwbin::eatWhiteSpace(parsed_line_feed[0]));
						std::string val(pwbin::eatWhiteSpace(parsed_line_feed[1]));
						// Found a group
						if(val.empty()){
								std::string group_name(parsed_line_feed[0]);
								group_name = pwbin::eatWhiteSpace(group_name);
								// Calculate amount of line whitespace to start string to
								// determine parent/child relationship
								level = pwbin::countFirstChar(line_feed," \t");
//								std::cout << "GROUP: " << group_name << ", LEVEL: " << level << std::endl;
								while(level <= levels.back() && parents.size() > 0){
										parents.pop_back();
										levels.pop_back();
								}
								levels.push_back(level);
								ParamBin* bin = parents.back();
								bin->setBin(group_name,std::unique_ptr<ParamBin>(new ParamBin));
								parents.push_back(bin->m_children[group_name].get());
						}
						else
                parents.back()->set(NamedParam<std::string>(name,val));
				} else {
						throw std::exception(); 
				}
    }
}

void ParamBin::loadParamFile(std::string FILE)
{
		namespace fs = std::filesystem;
		fs::path local_path = fs::path(FILE);
		if(!fs::exists(local_path)){
				fs::path current_dir = fs::current_path();
				fs::path full_path = current_dir / FILE;
				FILE = full_path.string();
		}
    std::ifstream fin(FILE);
    if(!fin.is_open())
    {
        fin.clear();
        throw ParamBinFileOpenException(FILE);
    }
		try{
		  	scanYAML(fin);
		} catch(...){
				fin.close();
        throw ParamBinFileReadException(FILE);
		}
}

void ParamBin::printBin(std::ostream& os) const{

    static const std::string EMPTY_CHARS = std::string(2,' ');
    static int depth = 0;
    os << std::endl;
    auto it = m_params.cbegin();
    while(it != m_params.cend()){
        std::string name("");
        for(int i = 0; i < depth; i++)
            name += EMPTY_CHARS;
        name += it->first;
        auto ait = m_ralias_map.find(it->first);
        if(ait != m_ralias_map.cend())
            name += '(' + ait->second + ')'; 

        os << std::setiosflags(std::ios::left) << std::setw(40) << name + ":";
        os << std::setw(16) << it->second << std::endl;
        it++;
    }
    os << std::endl;
    auto bit = m_children.cbegin();
    while(bit != m_children.cend()){
        std::string group_name("");
        for(int i = 0; i < depth; i++)
            group_name += EMPTY_CHARS;
        group_name += bit->first;
				group_name += ':';
        os << group_name << std::endl;
        depth++;
        bit->second->printBin(os);
        depth--;
        bit++;
    }
}

std::ostream& operator<<(std::ostream& os,const ParamBin& bin) 
{
    bin.printBin(os);
    return os;
}

bool ParamBin::inBin(const std::string& name) const
{
    if(m_params.count(name) > 0 || m_children.count(name) > 0) 
        return true;
    return false;
}

int ParamBin::size(const std::string& name) const
{
    auto it = m_params.find(name);
    if(it != m_params.cend())
        return pwbin::countCharacters(it->second,',') + 1;
    return 0;
}

bool ParamBin::empty() const
{
    return (m_params.empty() && m_children.empty() ? true : false);
}

ParamBin& ParamBin::getBin(const std::string& name) 
{
    auto it = m_children.find(name);
    if(it != m_children.cend())
        return *it->second;
    else
        throw ParamBinKeyException(name);
}

const ParamBin& ParamBin::getBin(const std::string& name) const
{
    auto it = m_children.find(name);
    if(it != m_children.cend())
        return *it->second;
    else
        throw ParamBinKeyException(name);
}


bool ParamBin::clear(const std::string& name)
{
    if(inBin(name)){
        m_params.erase(name);
        m_children.erase(name);
        return true;
    }
    else 
        return false;
}

void ParamBin::clear() 
{
    m_params.clear();
    m_children.clear();
    m_alias_map.clear();
    m_ralias_map.clear();
}

void ParamBin::setBool(const std::string& name,bool val){
    set(NamedParam<std::string>(name,(val ? "on" : "off")));
}

void ParamBin::set(const std::string& name,const ParamBin& bin)
{
    setBin(name,bin);
} 
 
// Pass by value requires constructing a new ParamBin with values copied in
void ParamBin::setBin(const std::string& name,const ParamBin& bin)
{
    std::unique_ptr<ParamBin> child_bin(new ParamBin(bin));
    child_bin->m_parent = this;
    m_children[name] = std::move(child_bin);
}

// Pass by pointers simply transfers ownership of pointer
void ParamBin::setBin(const std::string& name,std::unique_ptr<ParamBin> bin)
{
    bin->m_parent = this;
    m_children[name] = std::move(bin);
}

std::vector<std::string> ParamBin::inBin(const std::vector<std::string>& nameVec) const
{
    std::vector<std::string> good_names;
    for(unsigned int i = 0; i < nameVec.size(); i++){
        if(inBin(nameVec[i]))
            good_names.push_back(nameVec[i]);
    }
    return good_names;
}


std::vector<std::string> ParamBin::notInBin(const std::vector<std::string>& nameVec) const
{
    std::vector<std::string> bad_names;
    for(unsigned int i = 0; i < nameVec.size(); i++){
        if(!inBin(nameVec[i]))
            bad_names.push_back(nameVec[i]);
    }
    return bad_names;
}




