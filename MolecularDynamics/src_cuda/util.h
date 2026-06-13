#pragma once
#include <vector>
#include <map>

// std::map<std::string,double> readConfig(const std::string& filename){
//     std::ifstream iFile(filename);
//     std::map<std::string,double> config;
//     if (!iFile.is_open()){
//         throw std::runtime_error("Error: could not open config file - "+filename);
        
//     }
//     std::string line;
//     while (std::getline(iFile,line)){
//         std::stringstream ss(line);
//         std::string key;
//         std::string value;
//         if (std::getline(ss,key,'=') && std::getline(ss,value)){
//             config[key] = std::stod(value);
//         } 
//     }
//     return config;
// }

void pvdInit(std::ofstream& pvd){
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd << "  <Collection>\n";
}

void pvdWriteVTU(std::ofstream& pvd, int count, std::string fname){
    pvd<<"    <DataSet timestep=\"" << count << "\" file=\"" << fname << "\"/>\n";
}
