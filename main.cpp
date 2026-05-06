#include <iostream>
#include <string>
#include <fstream>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>

// JSON Interface module
#include "interface.h"

int main(int argc, char* argv[]) {
    // Check parameters, allowing an optional secondary argument for output filename
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_pdb_file> [output_file.json]\n";
        return 1;
    }

    std::string filepath = argv[1];
    std::string out_filepath;
    
    // Automatically generate output filename if not provided
    if (argc >= 3) {
        out_filepath = argv[2];
    } else {
        size_t last_slash = filepath.find_last_of("/\\");
        size_t last_dot = filepath.find_last_of(".");
        std::string base_name = filepath;
        
        if (last_slash != std::string::npos) {
            base_name = filepath.substr(last_slash + 1);
        }
        if (last_dot != std::string::npos && last_dot > last_slash) {
            base_name = base_name.substr(0, base_name.find_last_of("."));
        }
        out_filepath = base_name + "_xhpi.json";
    }

    try {
        std::cout << "[INFO] Loading structure: " << filepath << " ...\n";
        gemmi::Structure st = gemmi::read_structure_file(filepath);
        
        std::cout << "[INFO] Structure loaded. Scanning for interactions...\n";
        
        // Invoke the JSON interface to obtain the analysis string
        std::string json_output = xhpi::detect_xhpi_interactions_json(st);
        
        // Write the resulting JSON string to the output file
        std::ofstream outfile(out_filepath);
        if (!outfile.is_open()) {
            std::cerr << "[ERROR] Could not open output file: " << out_filepath << "\n";
            return 1;
        }
        
        outfile << json_output;
        outfile.close();
        
        std::cout << "[SUCCESS] Analysis complete. Output saved to: " << out_filepath << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "[FATAL] Exception occurred: " << e.what() << "\n";
        return 1;
    }

    return 0;
}