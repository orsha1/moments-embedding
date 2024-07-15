#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Define a structure to hold the data for each atom
struct Atom {
    int atomicNumber;
    float x;
    float y;
    float z;
};

// Function to read the CSV file
std::vector<Atom> readCSV(const std::string& filename) {
    std::vector<Atom> atoms;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return atoms;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        Atom atom;

        // Read the atomic number
        std::getline(ss, item, ',');
        atom.atomicNumber = std::stoi(item);

        // Read the x coordinate
        std::getline(ss, item, ',');
        atom.x = std::stof(item);

        // Read the y coordinate
        std::getline(ss, item, ',');
        atom.y = std::stof(item);

        // Read the z coordinate
        std::getline(ss, item, ',');
        atom.z = std::stof(item);

        atoms.push_back(atom);
    }

    file.close();
    return atoms;
}

int main() {
    std::string filename = "crystal_coordinates.csv";
    std::vector<Atom> atoms = readCSV(filename);

    // Print the atoms to verify the function works
    for (const Atom& atom : atoms) {
        std::cout << "Atomic Number: " << atom.atomicNumber
                  << ", x: " << atom.x
                  << ", y: " << atom.y
                  << ", z: " << atom.z << std::endl;
    }

    return 0;
}
