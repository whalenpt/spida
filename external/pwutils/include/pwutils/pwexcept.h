#ifndef PWEXCEPT_H
#define PWEXCEPT_H

#include <string>
#include <exception>

namespace pw{

class Exception : public std::exception
{
        public:
            explicit Exception(const std::string& message) : msg(message) {}
            explicit Exception(const std::string& where,const std::string& message)  
                {msg = "Error in " + where + ". " + message + ".";} 
            const char* what() const noexcept override {
                return msg.c_str();
            }
        private:
            std::string msg;
};


}

#endif



