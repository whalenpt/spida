
#ifndef PWDIR_H_
#define PWDIR_H_ 

#include <string>
#include <filesystem>

namespace pw{

class DirAux{

    public:
        DirAux();
        ~DirAux() {}
        std::string getCurrentDir() {return current_dir.string();}
        std::filesystem::path stringToPath(std::string dir_name);
        void createDirectory(std::filesystem::path dir_path,bool overwrite=true);
        void createDirectory(std::string dir_name,bool overwrite=true);
        bool checkDirectoryExists(std::filesystem::path dir_path);
        bool checkDirectoryExists(std::string dir_name);
        std::filesystem::path addTargetToDir(std::string target_name,std::filesystem::path dir_path);
        std::filesystem::path addTargetToDir(std::string target_name,std::string dir_name);
    private:
				std::filesystem::path current_dir;
};






}

#endif



