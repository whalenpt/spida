
#include "pwutils/report/json.hpp"
#include "pwutils/report/reporthelper.h"
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

namespace json{


void streamToJSON(std::ofstream& os,const pw::metadataMap& str_map) 
{
    pw::metadataMap::const_iterator it;
    for(it = str_map.begin(); it!= str_map.end(); it++){
        writeJSONLabel(os,(*it).first);
        writeJSONValue(os,(*it).second);
    }
}

void writeJSONLabel(std::ofstream& os,const std::string& label, const std::string& indent)
{
	os << indent << "\"" << label <<  "\" : ";
}

void writeJSONValue(std::ofstream& os,const std::string& value,
		const std::string& indent,bool end_value)
{
	std::string end_string = end_value ?  "\"" : "\","; 
    os << indent << "\"" << value << end_string << std::endl;
}

void mapToJSON(std::ofstream& os,const pw::metadataMap& str_map,bool end_value) 
{
	pw::metadataMap::const_iterator it;
	for(it = str_map.begin(); it!= str_map.end(); it++){
		writeJSONLabel(os,(*it).first);
		if(std::next(it) == str_map.end())
			writeJSONValue(os,(*it).second,"",end_value);
		else
			writeJSONValue(os,(*it).second,"",false);
	}
}




}







