
#ifndef PARAMBIN_HPP_ 
#define PARAMBIN_HPP_

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <exception>
#include <memory>
#include "ParamBin/parambinstrings.h"

std::ifstream& readNextLine(std::ifstream& fin,std::string& line_feed);
void lineToNameVal(const std::string& line_feed,std::string& name,std::string& vals);

class ParamBinException : public std::exception
{
    public:
        ParamBinException() {}
        virtual ~ParamBinException() {};
        virtual const char* what() const noexcept = 0;
};

class ParamBinKeyException : public ParamBinException
{
    public:
        explicit ParamBinKeyException(const std::string& nm) {
            msg = "The key " + nm + " was not found in the ParamBin.";
        }
        ~ParamBinKeyException() {};
        const char* what() const noexcept override {
            return msg.c_str();
        }
    private:
				std::string msg;
};

class ParamBinFileOpenException : public ParamBinException
{
    public:
        explicit ParamBinFileOpenException(const std::string& nm) {
            msg = "Failed to open the ParamBin parameter file " + nm + ".";
        }
        ~ParamBinFileOpenException() {};
        const char* what() const noexcept override {
            return msg.c_str();
        }
    private:
				std::string msg;
};

class ParamBinFileReadException : public ParamBinException
{
    public:
        explicit ParamBinFileReadException(const std::string& nm) {
            msg = "Failed to read the ParamBin parameter file " + nm + ".";
        }
        ~ParamBinFileReadException() {};
        const char* what() const noexcept override {
            return msg.c_str();
        }
    private:
				std::string msg;
};

class ParamBinScaleException : public ParamBinException
{
    public:
        explicit ParamBinScaleException(const std::string& scale){
            msg = "The scale " + scale + " could not be processed.";
        } 
        ~ParamBinScaleException() {};
        const char* what() const noexcept override {
            return msg.c_str();
        }
    private:
				std::string msg;
};


class ParamBin;
class NamedBin;
using ParamMap = std::map<std::string,std::string>; 
using AliasMap = std::map<std::string,std::string>;
using BinMap = std::map<std::string,std::unique_ptr<ParamBin>>;

template <class T>
class NamedParam
{
    public:
        NamedParam(const std::string& name,T val) :
            m_name(name), m_val(val) {}
        std::string name() const {return m_name;}
        T value() const {return m_val;}
    private:
        std::string m_name;
        T m_val;
};


class ParamBin{

    public:
        ParamBin();
        ~ParamBin() {}; 
        ParamBin(const ParamBin& bin);
        ParamBin(const char* FILE);
        ParamBin(std::string fileName);
        ParamBin& operator = (const ParamBin& bin);

        void loadParamFile(std::string FILE); 
        void loadParamFile(const char* FILE) {loadParamFile(std::string(FILE));}
  
        inline double getDbl(const std::string& nm) const { return get<double>(nm); }
        inline float getFloat(const std::string& nm) const { return get<float>(nm); }
        inline int getInt(const std::string& nm) const {return get<int>(nm);}
        inline std::string getStr(const std::string& nm) const {return get<std::string>(nm);}
        inline std::string getStrL(const std::string& nm) const {
            return pwbin::stringLowerCase(getStr(nm));}
        inline std::string getStrU(const std::string& nm) const {
            return pwbin::stringUpperCase(getStr(nm));}
        inline char getChar(const std::string& nm) const {return get<char>(nm);}
        inline bool getBool(const std::string& nm) const {
            std::string val = get<std::string>(nm);
            return (val == "on" ? true : false);
        }
        inline bool getBoolF(const std::string& nm) const {
            return (inBin(nm) ? getBool(nm) : false); 
        }
        inline bool getBoolT(const std::string& nm) const {
            return (inBin(nm) ? getBool(nm) : true); 
        }

        // get a child bin
        ParamBin& getBin(const std::string& name);
        const ParamBin& getBin(const std::string& name) const;

        // get the parent bin
        ParamBin& parentBin() {return *m_parent;} 
        const ParamBin& getParentBin() const {return *m_parent;}

        // get map of all params
        ParamMap getParamMap() const {return m_params;} 

        // get map of all child bins
        const BinMap& getChildBins() const {return m_children;}

        inline std::vector<double> getDblVec(const std::string& nm) const
        { return getVec<double>(nm); }
        inline std::vector<int> getIntVec(const std::string& nm) const
        { return getVec<int>(nm);}
        inline std::vector<float> getFloatVec(const std::string& nm) const
        { return getVec<float>(nm); }
        inline std::vector<std::string> getStrVec(const std::string& nm) const
        { return getVec<std::string>(nm); }
  
        template<typename T>
        void get(const std::string& nm,T& val) const;  // get a param of type T 

        template<typename T>
        T get(const std::string& nm) const; // get a param using get<T>(name)

        template<typename T>
        void get(const std::string& nm,std::vector<T>& val) const;  

        template<typename T>
        std::vector<T> getVec(const std::string& nm) const;  // get a vector using getVec<T>

        template<typename T>
        void set(const NamedParam<T>& named_param);

        template<typename T>
        void set(const std::string&,T val);

        template<typename T>
        void set(const std::string&,const std::vector<T>&);

        void set(const std::string&,const ParamBin& bin);

        void setBool(const std::string& nm,bool);

        void set(const NamedBin& named_bin);

        // Bins are deep copied into the tree structure
        void setBin(const std::string& nm,const ParamBin& bin);
        // Transfer ownership of bin 
        void setBin(const std::string& name,std::unique_ptr<ParamBin> bin);

        void setAlias(const std::string& name,const std::string& alias);

        friend std::ostream& operator<<(std::ostream&,const ParamBin& bin);

        template<typename T>
        ParamBin& operator<<(const NamedParam<T>& named_param);
        ParamBin& operator<<(const NamedBin& named_bin);

        bool inBin(const std::string&) const;
        std::vector<std::string> inBin(const std::vector<std::string>&) const;
        std::vector<std::string> notInBin(const std::vector<std::string>&) const;
        int size(const std::string& name) const;
        bool empty() const;
        void clear();
        bool clear(const std::string& name);

    private:

        ParamBin* m_parent;
        ParamMap m_params; 
        BinMap m_children;
        AliasMap m_alias_map;
        AliasMap m_ralias_map;

        std::string setParamKey(const std::string& name);
        std::string getStrParam(const std::string& name) const;
        void printBin(std::ostream& os) const;
        void scanYAML(std::ifstream& fin);

        bool searchAliasTree(const std::string& key,std::string& strval) const;
        bool searchParamMap(const std::string& key,std::string& strval) const;

        bool rootAliasSearch(const std::string& alias_key,std::string& strval) const;
};

class NamedBin
{
    public:
        NamedBin(const std::string& name,const ParamBin& bin) :
            m_name(name), m_bin(bin) {}
        std::string name() const {return m_name;}
        ParamBin bin() const {return m_bin;}
    private:
        std::string m_name;
        ParamBin m_bin;
};

template<typename T>
ParamBin& ParamBin::operator<<(const NamedParam<T>& named_param)
{
    set(named_param);
    return *this;
}

template<typename T> 
std::string convertToString(T val) 
{
    // Use std::to_string approach to get string representation of value
	std::string str(std::to_string(val));
	// check for decimal number, else use built-in std::to_string
    if(str.find('.') == std::string::npos)
        return str;
   	int offset = 1;
	if(str.find_last_not_of('0')==str.find('.')){
		offset = 0;
	}
	str.erase(str.find_last_not_of('0') + offset,std::string::npos);


    // check for zeros (-> std::to_string does not provide high enough precision)
    bool to_string_worked = false;
    for(auto& ch : str){
        if(ch != '0' && ch != '.' && ch !='-'){
            to_string_worked = true;
            break;
        }
    }

    // Use stringstream approach to get string representation of value
	std::ostringstream oss;
	oss.precision(16);
	oss << std::scientific << val;
	std::string stre = oss.str();

    // std::to_string did not work on decimal/scientific number -> must use stringstream
	if(!to_string_worked)
	    str = stre;


    std::string decimal_part = pwbin::subString(stre,'.','e');
    size_t pos = decimal_part.find_last_not_of('0');
    if(pos != std::string::npos){
        decimal_part = decimal_part.substr(0,pos+1);
    } else{
            // 1.00000e5 - > 1.0e5
        decimal_part = "0";
    }
    std::string exponent_part = pwbin::stripLast(stre,'e');
    std::string exp_num,exp_sign;
    // Remove '+' exponent sign for string storage
    if(exponent_part.find_first_of('-') != std::string::npos){
        exp_sign = "-";
        exp_num = exponent_part.substr(1,std::string::npos);
    } else if(exponent_part.find_first_of('+') != std::string::npos){
        exp_sign = "";
        exp_num = exponent_part.substr(1,std::string::npos);
    } else{
        exp_sign = "";
        exp_num = exponent_part.substr(0,std::string::npos);
    }
    pos = exp_num.find_first_not_of('0');
    if(pos != std::string::npos){
    exp_num = exp_num.substr(pos,std::string::npos);
    } else
        // replace 1.3e000 with 1.3e0
        exp_num = "0";

    exponent_part = exp_sign + exp_num;

    std::string whole_part = pwbin::stripFirst(stre,'.');
    stre =  whole_part + "." + decimal_part + "e" + exponent_part;
    // std::cout << "EXPONENTIAL FORM: " << stre << std::endl;
	// Check if exponential form is significantly shorter than fixed point
	int FEWER_SPACES = 0;
	if( (stre.size()+FEWER_SPACES) < str.size()) 
        str = stre;
    return str;
}

template<typename T> 
std::string convertToString(const std::vector<T>& val) 
{
    if(val.empty())
        return "";

    std::string vecstr = convertToString(val[0]);
    for(unsigned int i = 1; i < val.size(); i++){
        std::string valstr = convertToString(val[i]);
        vecstr = vecstr + "," + valstr;
    }
    return vecstr;
}

template<> 
std::string convertToString(std::string val);

template<> 
std::string convertToString(const char* val);

template<> 
std::string convertToString(char val);

template <typename T>
void convertFromString(const std::string& str,T& val)
{
    std::stringstream ss;
    ss.setf(std::ios_base::showpoint);
    ss.precision(12);
    ss << str;
    ss >> val;
}

template <>
void convertFromString(const std::string& str,char& val);

template <>
void convertFromString(const std::string& str,std::string& val);

// ADD A PARAM TO THE PARAM BIN
template<typename T>
void ParamBin::set(const std::string& name,T val) 
{
    set(NamedParam<T>(name,val));
}

// ADD A PARAM SET TO THE PARAM BIN
template<typename T>
void ParamBin::set(const std::string& name,const std::vector<T>& val) 
{
    set(NamedParam<std::vector<T>>(name,val));
}

template<typename T>
void ParamBin::set(const NamedParam<T>& named_param)
{
    std::string key = setParamKey(named_param.name());
    std::string strVal = convertToString(named_param.value());
    m_params[key] = strVal;
}

template<typename T>
void ParamBin::get(const std::string& name,T& val) const
{
    std::string strval = getStrParam(name);
    convertFromString<T>(strval,val);
}

template<typename T>
T ParamBin::get(const std::string& name) const
{
    T val;
    get(name,val);
    return val;
}

template<typename T>
void ParamBin::get(const std::string& name,std::vector<T>& vals) const
{
    std::string strval = getStrParam(name);
    std::vector<std::string> strvec = pwbin::parseString(strval,',');
    vals.clear();
    for(auto strval : strvec){
        T val;
        convertFromString<T>(strval,val);
        vals.push_back(val);
    }
}

template<typename T>
std::vector<T> ParamBin::getVec(const std::string& name) const
{
    std::vector<T> vals;
    get(name,vals);
    return vals;
}



#endif



    
 

