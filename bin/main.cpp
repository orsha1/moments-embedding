#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <set>
#include <numeric>
#include <cmath>

struct Settings
{
    int n_pbc;
    int l_path;
    int l_neis;

    std::string datapath;
    std::string type;
    std::string runpath;
    std::string input_coordinates = "input_coordinates.csv";
    std::string input_lattice_parameters = "input_lattice_parameters.csv";
    std::string name = "step";

    bool embed;
    bool quick;
    double bond_thr;
    bool homoatomic_bonding;
};

struct Atom
{
    int symbol;
    double x;
    double y;
    double z;
    std::vector<std::vector<int>> self_returning_paths;
    std::vector<std::vector<int>> cycles;
    std::vector<std::vector<int>> neis;
};

class Structure
{
public:
    int nat;
    std::vector<Atom> atoms;
    std::vector<std::vector<double>> lattice_parameters;
    std::vector<std::vector<double>> coordinates;
    std::vector<std::vector<double>> relative_coordinates;
    std::vector<int> symbols;
    std::vector<std::vector<int>> image;
    std::vector<std::vector<int>> cycles_full;
    std::vector<std::vector<int>> self_returning_full;
    std::vector<std::vector<int>> neis_full;

    std::vector<std::vector<int>> cycles;
    std::vector<std::vector<int>> self_returning;
    std::vector<std::vector<int>> neis;

    std::vector<std::vector<int>> neighbors;
    std::vector<std::vector<double>> distances;

    std::unordered_map<int, std::vector<int>> cycles_bonds;
    std::unordered_map<int, std::vector<int>> paths_bonds;

    std::unordered_map<int, std::vector<double>> cycles_coords_step;
    std::unordered_map<int, std::vector<double>> paths_coords_step;

    std::unordered_map<int, double> cycles_distances_step;
    std::unordered_map<int, double> paths_distances_step;
};

bool vectorCompare(const std::vector<int> &a, const std::vector<int> &b)
/**
 * Compare two vectors for equality.
 *
 * @param a The first vector.
 * @param b The second vector.
 * @return True if vectors a and b are equal, false otherwise.
 */
{
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
}

std::vector<int> uniqueValues(const std::vector<int> &inputVector)
{
    // Create a copy of the input vector to avoid modifying the original
    std::vector<int> resultVector = inputVector;

    // Sort the vector to bring duplicate values together
    std::sort(resultVector.begin(), resultVector.end());

    // Use unique algorithm to remove adjacent duplicates
    auto uniqueEnd = std::unique(resultVector.begin(), resultVector.end());

    // Resize the vector to keep only the unique values
    resultVector.resize(std::distance(resultVector.begin(), uniqueEnd));

    return resultVector;
}

void removeElements(std::vector<int> mainVector, std::vector<int> elementsToRemove)
{

    // Sort both vectors to enable efficient removal of common elements
    std::sort(mainVector.begin(), mainVector.end());
    std::sort(elementsToRemove.begin(), elementsToRemove.end());

    // Use std::set_difference to find elements to be removed
    auto it = std::set_difference(mainVector.begin(), mainVector.end(),
                                  elementsToRemove.begin(), elementsToRemove.end(),
                                  mainVector.begin());

    // Resize the vector to keep only the unique values
    mainVector.resize(std::distance(mainVector.begin(), it));
}

Settings loadSettings(const std::string &filename)
{
    /**
     * Load settings from a configuration file.
     *
     * @param filename The name of the configuration file.
     * @return An instance of the Settings struct populated with values from the configuration file.
     */
    Settings settings; // Create an instance of the Settings struct

    std::ifstream file(filename);

    // Check if the file is open
    if (file.is_open())
    {
        std::string line;

        // Read the file line by line
        while (std::getline(file, line))
        {
            // Skip comments and empty lines
            if (!line.empty() && line[0] != '#')
            {
                std::istringstream iss(line);
                std::string key, value;

                // Split each line into key and value
                if (std::getline(iss, key, '=') && std::getline(iss, value))
                {
                    // Check the key and set the corresponding member in the struct
                    if (key == "datapath")
                    {
                        settings.datapath = value;
                    }
                    else if (key == "n_pbc")
                    {
                        settings.n_pbc = std::stoi(value);
                    }
                    else if (key == "bond_thr")
                    {
                        settings.bond_thr = std::stod(value);
                    }
                    else if (key == "l_path")
                    {
                        settings.l_path = std::stoi(value);
                    }
                    else if (key == "l_neis")
                    {
                        settings.l_neis = std::stoi(value);
                    }
                    else if (key == "embed")
                    {
                        std::istringstream(value) >> std::boolalpha >> settings.embed;
                    }
                    else if (key == "runpath")
                    {
                        settings.runpath = value;
                    }
                    else if (key == "type")
                    {
                        settings.type = value;
                    }
                    else if (key == "coordinates")
                    {
                        settings.input_coordinates = value;
                    }
                    else if (key == "lattice_parameters")
                    {
                        settings.input_lattice_parameters = value;
                    }
                    else if (key == "name")
                    {
                        settings.name = value;
                    }
                    else if (key == "quick")
                    {
                        std::istringstream(value) >> std::boolalpha >> settings.quick;
                    }
                    else if (key == "homoatomic_bonding")
                    {
                        std::istringstream(value) >> std::boolalpha >> settings.homoatomic_bonding;
                    }
                }
            }
        }

        file.close(); // Close the file after reading
    }
    else
    {
        // Print an error message if the file cannot be opened
        std::cerr << "Error opening file: " << filename << std::endl;
    }

    return settings; // Return the populated Settings struct
}

std::vector<Atom> readCoordsCSV(const std::string &filename)
{
    std::vector<Atom> data;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening CSV file: " << filename << std::endl;
        return {};
    }

    std::string line;

    // Read the file line by line
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        Atom row;

        // Read symbol, x, y, z from the line
        ss >> row.symbol;
        ss.ignore(); // Ignore the comma
        ss >> row.x;
        ss.ignore(); // Ignore the comma
        ss >> row.y;
        ss.ignore(); // Ignore the comma
        ss >> row.z;

        // Add the row to the data vector
        data.push_back(row);
    }

    file.close();

    return data;
}

std::vector<std::vector<double>> readLatticeParamsCSV(const std::string &filename)
{
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening CSV file: " << filename << std::endl;
        return {};
    }

    std::string line;

    // Read the file line by line
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::vector<double> row;

        double value;
        // Read symbol, x, y, z from the line
        while (ss >> value)
        {
            row.push_back(value);
            if (ss.peek() == ',')
                ss.ignore();
        }

        // Add the row to the data vector
        data.push_back(row);
    }

    file.close();

    return data;
}

std::vector<std::vector<int>> readPaths(const std::string &filename)
{
    std::ifstream file(filename);
    std::vector<std::vector<int>> data;

    if (file.is_open())
    {
        std::string line;

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::vector<int> row;

            int value;
            while (iss >> value)
            {
                row.push_back(value);
            }

            data.push_back(row);
        }

        file.close();
    }
    else
    {
        std::cerr << "Error opening file: " << filename << std::endl;
    }

    return data;
}

std::unordered_map<int, std::vector<int>> readBonds(const std::string &filename)
{
    std::unordered_map<int, std::vector<int>> dataMap;

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return dataMap; // Return an empty map if the file cannot be opened
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        int key;
        iss >> key;

        std::vector<int> values;
        int value;
        while (iss >> value)
        {
            values.push_back(value);
        }

        dataMap[key] = values;
    }

    file.close();
    return dataMap;
}

std::vector<std::vector<int>> readBinaryFile(const std::string &filename)
{
    std::ifstream file(filename, std::ios::binary | std::ios::in);

    if (!file.is_open())
    {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return {};
    }

    // Read the number of rows
    std::size_t rows;
    file.read(reinterpret_cast<char *>(&rows), sizeof(std::size_t));

    // Read each row
    std::vector<std::vector<int>> data(rows);
    for (auto &row : data)
    {
        // Read the size of the inner vector
        std::size_t innerSize;
        file.read(reinterpret_cast<char *>(&innerSize), sizeof(std::size_t));
        // Resize the inner vector and read its elements
        row.resize(innerSize);
        file.read(reinterpret_cast<char *>(row.data()), innerSize * sizeof(int));
    }

    file.close();
    return data;
}

/**
 * Multiply two matrices.
 *
 * @param mat1 The first matrix.
 * @param mat2 The second matrix.
 * @return The result of multiplying mat1 and mat2.
 *         If the matrices cannot be multiplied, an empty vector is returned.
 */
std::vector<std::vector<double>> matrixMultiply(const std::vector<std::vector<double>> &mat1, const std::vector<std::vector<double>> &mat2)
{
    size_t rows1 = mat1.size();
    size_t cols1 = mat1[0].size();
    size_t rows2 = mat2.size();
    size_t cols2 = mat2[0].size();

    if (cols1 != rows2)
    {
        std::cerr << "Error: Matrices cannot be multiplied. Invalid dimensions." << std::endl;
        return std::vector<std::vector<double>>();
    }

    std::vector<std::vector<double>> result(rows1, std::vector<double>(cols2, 0));

    for (size_t i = 0; i < rows1; ++i)
    {
        for (size_t j = 0; j < cols2; ++j)
        {
            for (size_t k = 0; k < cols1; ++k)
            {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return result;
}

/**
 * Calculate the pairwise distances between two sets of coordinates.
 *
 * @param coordinates_i The coordinates of the first set.
 * @param coordinates_j The coordinates of the second set.
 * @return A matrix containing pairwise distances between coordinates_i and coordinates_j.
 */
std::vector<std::vector<double>> cdist(std::vector<std::vector<double>> &coordinates_i, std::vector<std::vector<double>> &coordinates_j)
{
    std::vector<std::vector<double>> distances(coordinates_i.size(), std::vector<double>(coordinates_j.size()));

    for (int i = 0; i < coordinates_i.size(); ++i)
    {
        double x_i = coordinates_i[i][0];
        double y_i = coordinates_i[i][1];
        double z_i = coordinates_i[i][2];

        for (int j = 0; j < coordinates_j.size(); ++j)
        {
            double dx = x_i - coordinates_j[j][0];
            double dy = y_i - coordinates_j[j][1];
            double dz = z_i - coordinates_j[j][2];

            distances[i][j] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
    return distances;
}

/**
 * Check if all elements in a vector are unique.
 *
 * @param vec The input vector.
 * @return True if all elements are unique, false otherwise.
 */
bool areAllUnique(const std::vector<int> &vec)
{
    std::vector<int> tempVec = vec;

    std::sort(tempVec.begin(), tempVec.end());

    auto last = std::unique(tempVec.begin(), tempVec.end());

    return std::distance(tempVec.begin(), last) == tempVec.size();
}

// Function to write atomic structure data to XYZ and VASP files
void writeStructure(const Structure &print_structure, Settings settings, std::string filename)
{
    std::map<int, int> counts;

    std::ofstream xyzfile;
    xyzfile.open(settings.datapath + filename + ".xyz");

    // Write atomic coordinates to XYZ file
    for (int row = 0; row < print_structure.nat; ++row)
    {
        xyzfile << print_structure.symbols[row] << " ";
        counts[print_structure.symbols[row]]++;

        for (int col = 0; col < 3; ++col)
        {
            xyzfile << std::fixed << std::setprecision(8) << print_structure.coordinates[row][col] << " ";
        }
        xyzfile << std::endl;
    }

    // Map chemical symbols to atomic numbers
    std::map<int, std::string> chemical_symbols = {
        {42, "Mo"},
        {56, "Ba"},
        {22, "Ti"},
        {8, "O"}};

    std::ofstream vaspfile;
    vaspfile.open(settings.datapath + filename + ".vasp");
    vaspfile << "New Structure" << std::endl;
    vaspfile << 1 << std::endl;

    // Write lattice parameters to VASP file
    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            vaspfile << std::fixed << std::setprecision(8) << print_structure.lattice_parameters[row][col] << " ";
        }
        vaspfile << std::endl;
    }

    // Write chemical symbols to VASP file
    for (const auto &pair : counts)
    {
        vaspfile << chemical_symbols[pair.first] << " ";
    }
    vaspfile << std::endl;

    // Write counts of each chemical symbol to VASP file
    for (const auto &pair : counts)
    {
        vaspfile << pair.second << " ";
    }
    vaspfile << std::endl;
    vaspfile << "Cartesian" << std::endl;

    // Write atomic coordinates (Cartesian) to VASP file
    for (const auto &pair : counts)
    {
        for (int row = 0; row < print_structure.nat; ++row)
        {
            if (print_structure.symbols[row] == pair.first)
            {
                for (int col = 0; col < 3; ++col)
                {
                    vaspfile << std::fixed << std::setprecision(8) << print_structure.coordinates[row][col] << " ";
                }
                vaspfile << std::endl;
            }
        }
    }
}

// Function to write embedding data to output files
void writeEmbeddingData(const Structure &structure, const Settings &settings)
{
    std::ofstream neighbors_out;
    std::ofstream neighbors_out_translated;

    std::ofstream cycles_out;
    std::ofstream cycles_out_translated;

    std::ofstream self_returning_paths_out;
    std::ofstream self_returning_paths_out_translated;

    std::ofstream neis_out;
    std::ofstream neis_out_translated;

    neighbors_out.open(settings.datapath + "neighbors_pbc.dat");
    neighbors_out_translated.open(settings.datapath + "neighbors_pbc_translated.dat");

    cycles_out.open(settings.datapath + "cycles_pbc.dat");
    cycles_out_translated.open(settings.datapath + "cycles_pbc_translated.dat");

    self_returning_paths_out.open(settings.datapath + "self_returning_pbc.dat");
    self_returning_paths_out_translated.open(settings.datapath + "self_returning_pbc_translated.dat");

    neis_out.open(settings.datapath + "neis.dat");
    neis_out_translated.open(settings.datapath + "neis_translated.dat");

    for (int i = 0; i < structure.nat; ++i)
    {
        std::vector<int> image = structure.image[i];

        // Write neighbor information to files
        neighbors_out << i << " " << structure.symbols[i] << " " << structure.atoms[i].x << " " << structure.atoms[i].y << " " << structure.atoms[i].z << " ";
        neighbors_out_translated << i << " " << structure.image[i][0] << " " << structure.image[i][1] << " " << structure.image[i][2] << " " << structure.image[i][3] << " ";

        for (int j = 0; j < structure.neighbors[i].size(); ++j)
        {
            neighbors_out << structure.neighbors[i][j] << " ";
            for (int k = 0; k < 4; ++k)
            {
                neighbors_out_translated << structure.image[structure.neighbors[i][j]][k] << " ";
            }
        }
        neighbors_out << std::endl;
        neighbors_out_translated << std::endl;

        // Write cycle information to files
        if ((image[0] == 0) && (image[1] == 0) && (image[2] == 0))
        {
            for (int j = 0; j < structure.atoms[i].cycles.size(); ++j)
            {
                cycles_out_translated << structure.image[i][3] << " " << structure.symbols[i] << " " << structure.atoms[i].cycles[j].size() << " ";
                cycles_out << i << " " << structure.symbols[i] << " " << structure.atoms[i].cycles[j].size() << " ";

                for (int m = 0; m < structure.atoms[i].cycles[j].size(); ++m)
                {
                    cycles_out << structure.atoms[i].cycles[j][m] << " ";
                    for (int k = 0; k < 4; ++k)
                    {
                        cycles_out_translated << structure.image[structure.atoms[i].cycles[j][m]][k] << " ";
                    }
                }
                cycles_out << std::endl;
                cycles_out_translated << std::endl;
            }
        }
        // Write neis information to files
        if ((image[0] == 0) && (image[1] == 0) && (image[2] == 0))
        {
            for (int j = 0; j < structure.atoms[i].neis.size(); ++j)
            {
                neis_out_translated << structure.image[i][3] << " " << structure.symbols[i] << " " << j << " ";
                neis_out << i << " " << structure.symbols[i] << " " << structure.atoms[i].neis[j].size() << " ";

                for (int m = 0; m < structure.atoms[i].neis[j].size(); ++m)
                {
                    neis_out << structure.atoms[i].neis[j][m] << " ";
                    for (int k = 0; k < 4; ++k)
                    {
                        neis_out_translated << structure.image[structure.atoms[i].neis[j][m]][k] << " ";
                    }
                }
                neis_out << std::endl;
                neis_out_translated << std::endl;
            }
        }
        // Write self-returning path information to files
        if ((image[0] == 0) && (image[1] == 0) && (image[2] == 0))
        {
            for (int j = 0; j < structure.atoms[i].self_returning_paths.size(); ++j)
            {
                self_returning_paths_out_translated << structure.image[i][3] << " " << structure.symbols[i] << " " << structure.atoms[i].self_returning_paths[j].size() << " ";
                self_returning_paths_out << i << " " << structure.symbols[i] << " " << structure.atoms[i].self_returning_paths[j].size() << " ";

                for (int m = 0; m < structure.atoms[i].self_returning_paths[j].size(); ++m)
                {
                    self_returning_paths_out << structure.atoms[i].self_returning_paths[j][m] << " ";
                    for (int k = 0; k < 4; ++k)
                    {
                        self_returning_paths_out_translated << structure.image[structure.atoms[i].self_returning_paths[j][m]][k] << " ";
                    }
                }
                self_returning_paths_out << std::endl;
                self_returning_paths_out_translated << std::endl;
            }
        }
    }
}

// Function to write bonding data to output files
void writeBondingData(const Structure &structure, const Settings &settings)
{
    std::ofstream bonds_cycles_out_translated;
    std::ofstream bonds_self_returning_paths_out_translated;

    bonds_cycles_out_translated.open(settings.datapath + "bonds_cycles_out_translated.dat");
    bonds_self_returning_paths_out_translated.open(settings.datapath + "bonds_self_returning_paths_out_translated.dat");

    // Write cycle bonds to translated file
    for (const auto &entry : structure.cycles_bonds)
    {
        int key = entry.first;
        const std::vector<int> &values = entry.second;
        bonds_self_returning_paths_out_translated << key << " ";
        for (int value : values)
        {
            bonds_self_returning_paths_out_translated << value << " ";
        }
        bonds_self_returning_paths_out_translated << std::endl;
    }

    // Write self-returning path bonds to translated file
    for (const auto &entry : structure.paths_bonds)
    {
        int key = entry.first;
        const std::vector<int> &values = entry.second;
        bonds_cycles_out_translated << key << " ";
        for (int value : values)
        {
            bonds_cycles_out_translated << value << " ";
        }
        bonds_cycles_out_translated << std::endl;
    }
}

void writeBinaryFile(const std::string &filename, const std::vector<std::vector<int>> &data)
{
    std::ofstream file(filename, std::ios::binary | std::ios::out);

    if (!file.is_open())
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    // Write the number of rows
    std::size_t rows = data.size();
    file.write(reinterpret_cast<char *>(&rows), sizeof(std::size_t));

    // Write each row
    for (const auto &row : data)
    {
        // Write the size of the inner vector
        std::size_t innerSize = row.size();
        file.write(reinterpret_cast<char *>(&innerSize), sizeof(std::size_t));
        // Write the elements of the inner vector
        file.write(reinterpret_cast<const char *>(row.data()), innerSize * sizeof(int));
    }

    file.close();
}

std::vector<double> pbc_d(const std::vector<std::vector<double>> &lp, const std::vector<int> &v)
{
    {

        int n = lp[0].size(); // Assuming all rows of lp have the same size
        std::vector<double> d(n, 0.0);

        for (int i = 0; i < 3; ++i)
        { // Assuming lp always has 3 rows
            for (int j = 0; j < n; ++j)
            {
                d[j] += lp[i][j] * v[i];
            }
        }

        return d;
    }
}

Structure PBC(const Structure &original_structure, const int &n)
{
    /**
     * Periodic Boundary Condition (PBC) Expansion.
     *
     * @param original_structure The original structure.
     * @param n The expansion factor.
     * @return A new structure with coordinates expanded based on PBC.
     */
    Structure structure_pbc;
    std::vector<std::vector<int>> image;
    for (int i = -n; i <= n; ++i)
    {
        for (int j = -n; j <= n; ++j)
        {
            for (int k = -n; k <= n; ++k)
            {
                // std::vector<std::vector<double>> displace_vector = original_structure.lattice_parameters;
                std::vector<int> image_t = {i, j, k};
                std::vector<double> displace_vector = pbc_d(original_structure.lattice_parameters, image_t);
                std::vector<std::vector<double>> displaced_coordinates = original_structure.coordinates;
                // std::cout << "Displace vector " << displace_vector[0] << " " << displace_vector[1] << " " << displace_vector[2] << std::endl;
                for (int row = 0; row < original_structure.nat; ++row)
                {
                    // std::cout << "Coordinates " << displaced_coordinates[row][0] << " " << displaced_coordinates[row][1] << " " << displaced_coordinates[row][2] << std::endl;
                    image.push_back({i, j, k, row});
                    for (int col = 0; col < 3; ++col)
                    {
                        displaced_coordinates[row][col] = displaced_coordinates[row][col] + displace_vector[col];
                    }
                    // std::cout << "Displaced Coordinates " << displaced_coordinates[row][0] << " " << displaced_coordinates[row][1] << " " << displaced_coordinates[row][2] << std::endl;
                    structure_pbc.coordinates.push_back(displaced_coordinates[row]);
                    Atom atom_i = {original_structure.symbols[row], displaced_coordinates[row][0], displaced_coordinates[row][1], displaced_coordinates[row][2]};
                    structure_pbc.atoms.push_back(atom_i);
                    structure_pbc.symbols.push_back(original_structure.symbols[row]);
                }
            }
        }
    }
    // -1.86, 4.2, 5.76
    // -16.20, -9.663, -15.642
    // 7.4368090630,0.0000000000,0.0000000000
    // 0.0000000000,7.7957539558,0.0000000000
    // -7.4039145041,0.0000000000,8.0082301811
    std::vector<std::vector<double>> lattice_parameters = original_structure.lattice_parameters;

    for (size_t i = 0; i < lattice_parameters.size(); ++i)
    {
        for (size_t j = 0; j < lattice_parameters[i].size(); ++j)
        {
            lattice_parameters[i][j] = (2 * n + 1) * lattice_parameters[i][j];
        }
    }
    structure_pbc.lattice_parameters = lattice_parameters;
    structure_pbc.nat = structure_pbc.coordinates.size();
    structure_pbc.image = image;
    return structure_pbc;
}

std::vector<std::vector<int>> mapNeighbors(std::vector<std::vector<double>> &distances, std::vector<int> &symbols, double thr, bool homoatomic_bonding)
{
    /**
     * Map neighboring atoms based on a distance threshold.
     *
     * @param distances Pairwise distances between atoms.
     * @param thr The distance threshold.
     * @return A vector of vectors representing neighbors for each atom.
     */
    std::vector<std::vector<int>> neighbors;

    for (int i = 0; i < distances.size(); ++i)
    {
        std::vector<int> neighbors_i;
        if (homoatomic_bonding)
        {
            // Homoatomic bonding: consider bonds between atoms of the same type
            for (int j = 0; j < distances[i].size(); ++j)
            {

                if (distances[i][j] < thr && distances[i][j] > 0)
                {

                    neighbors_i.push_back(j);
                }
            }
        }
        else
        {
            // Heteroatomic bonding only: consider bonds between atoms of different types
            int symbol_i = symbols[i];
            for (int j = 0; j < distances[i].size(); ++j)
            {
                int symbol_j = symbols[j];
                if (distances[i][j] < thr && distances[i][j] > 0 && symbol_i != symbol_j)
                {

                    neighbors_i.push_back(j);
                }
            }
        }
        neighbors.push_back(std::move(neighbors_i));
    }
    return neighbors;
}

std::vector<std::vector<int>> expandPaths(std::vector<std::vector<int>> neighbors_pbc, std::vector<std::vector<int>> paths)
{
    /**
     * Expand paths by appending neighboring atoms.
     *
     * @param neighbors_pbc Neighbors of atoms based on PBC.
     * @param paths Paths to be expanded.
     * @return New paths generated by appending neighboring atoms.
     */

    std::vector<std::vector<int>> new_paths;

    for (int i = 0; i < paths.size(); ++i)
    {
        std::vector<int> current_path = paths[i];
        int last_member = current_path[current_path.size() - 1];
        std::vector<int> last_member_neis = neighbors_pbc[last_member];

        for (int j = 0; j < last_member_neis.size(); ++j)
        {
            std::vector<int> path_ij = current_path;
            path_ij.push_back(last_member_neis[j]);
            new_paths.push_back(path_ij);
        }
    }
    return new_paths;
}

std::vector<std::vector<std::vector<int>>> generatePaths(std::vector<std::vector<int>> neighbors_pbc, int starting_point, int max_length, int max_length_neis)
{
    /**
     * Generate self-returning paths and cycles up to a specified length.
     *
     * @param neighbors_pbc Neighbors of atoms based on PBC.
     * @param starting_point The starting atom for path generation.
     * @param max_length The maximum length of paths to be generated.
     * @return A vector containing two types of paths: self-returning paths and cycles.
     */
    std::map<int, std::vector<std::vector<int>>> paths;
    paths[0] = {{starting_point}};
    std::vector<std::vector<int>> self_returning;
    std::vector<std::vector<int>> cycles;
    std::vector<std::vector<int>> neis;

    for (int i = 1; i <= max_length; ++i)
    {
        paths[i] = expandPaths(neighbors_pbc, paths[i - 1]);
    }

    for (int i = 1; i <= max_length_neis; ++i)
    {
        std::vector<int> neis_n;

        for (int j = 0; j < paths[i].size(); ++j)
        {
            std::vector<int> &path = paths[i][j];
            neis_n.push_back(path[path.size() - 1]);
        }
        std::vector<int> neis_n_clean = uniqueValues(neis_n);
        if (i > 0)
        {
            for (int k = 1; k < i; ++k)
            {
                removeElements(neis_n_clean, neis[k - 1]);
            }
        }
        neis.push_back(neis_n_clean);
    }

    for (int i = 0; i < paths.size(); ++i)
    {
        for (int j = 0; j < paths[i].size(); ++j)
        {
            std::vector<int> &path = paths[i][j];
            if (path[0] == path[path.size() - 1])
            {
                self_returning.push_back(path);
                if (areAllUnique(std::vector<int>(path.begin() + 1, path.end())))
                {
                    cycles.push_back(path);
                }
            }
        }
    }

    return std::vector<std::vector<std::vector<int>>>{self_returning, cycles, neis};
}

/**
 * Obtain bonds from translated paths.
 *
 * @param paths_translated Translated paths containing information about bonds.
 * @return An unordered map where keys are bond indices and values are pairs of atoms forming the bond.
 */
std::unordered_map<int, std::vector<int>> obtainBonds(std::vector<std::vector<int>> &paths_translated)
{
    std::unordered_map<int, std::vector<int>> bonds_map;
    std::vector<std::vector<int>> paths_translated_indexed;

    for (int i = 0; i < paths_translated.size(); ++i)
    {
        std::vector<int> cycle = paths_translated[i];
        std::vector<int> only_cyc(cycle.begin() + 3, cycle.end());

        int length = cycle[2];
        if (length == 1)
        {
            paths_translated_indexed.push_back(cycle);
        }
        else
        {
            std::vector<int> cycle_indexed(cycle.begin(), cycle.begin() + 3);

            for (int j = 0; j < length - 1; ++j)
            {
                std::vector<int> pair(only_cyc.begin() + j * 4, only_cyc.begin() + j * 4 + 8);

                int idx0 = 0;
                int idx1 = 0;
                for (int k = 0; k < 4; ++k)
                {
                    idx0 += (10 + pair[k]) * std::pow(10, k);
                    idx1 += (10 + pair[k + 4]) * std::pow(10, k);
                }

                cycle_indexed.push_back(idx0 + idx1);

                bonds_map[idx0 + idx1] = pair;
            }
            paths_translated_indexed.push_back(cycle_indexed);
        }
    }

    return bonds_map;
}

// Function to initialize the embedding process
std::vector<Structure> initEmbed(std::string inputname0, std::string inputname1, Settings settings)
{
    Structure structure;

    // Read atomic coordinates from CSV file
    structure.atoms = readCoordsCSV(inputname0);
    // Read lattice parameters from CSV file
    structure.lattice_parameters = readLatticeParamsCSV(inputname1);

    structure.nat = structure.atoms.size();

    // Populate relative coordinates and symbols
    for (int i = 0; i < structure.nat; ++i)
    {
        std::vector row = {structure.atoms[i].x,
                           structure.atoms[i].y,
                           structure.atoms[i].z};
        structure.relative_coordinates.push_back(row);
        structure.symbols.push_back(structure.atoms[i].symbol);
    }

    // Calculate absolute coordinates using matrix multiplication
    structure.coordinates = matrixMultiply(structure.relative_coordinates, structure.lattice_parameters);

    // Apply periodic boundary conditions
    Structure structure_pbc = PBC(structure, settings.n_pbc);

    // Write structures to files
    writeStructure(structure, settings, "structure");
    writeStructure(structure_pbc, settings, "structure_pbc");

    // Calculate distances and map neighbors for the PBC structure
    structure_pbc.distances = cdist(structure_pbc.coordinates, structure_pbc.coordinates);

    // for (int i = 0; i < structure_pbc.distances[1984].size(); ++i)
    //{
    //     if ((structure_pbc.image[i][0] == 0) & (structure_pbc.image[i][1] == 0) & (structure_pbc.image[i][2] == 0))
    //     {
    //         std::cout << i << " " << structure_pbc.symbols[i] << " " << structure_pbc.image[i][0] << " " << structure_pbc.image[i][1] << " " << structure_pbc.image[i][2] << " " << structure_pbc.coordinates[i][0] << " " << structure_pbc.coordinates[i][1] << " " << structure_pbc.coordinates[i][2] << " " << structure_pbc.distances[1984][i] << std::endl;
    //     }
    //     //if ((structure_pbc.image[i][0] == 1) & (structure_pbc.image[i][1] == 0) & (structure_pbc.image[i][2] == 0))
    //     //{
    //     //    std::cout << i << " " << structure_pbc.symbols[i] << " " << structure_pbc.image[i][0] << " " << structure_pbc.image[i][1] << " " << structure_pbc.image[i][2] << " " << structure_pbc.coordinates[i][0] << " " << structure_pbc.coordinates[i][1] << " " << structure_pbc.coordinates[i][2] << " " << structure_pbc.distances[1984][i] << std::endl;
    //     //}
    // }

    structure_pbc.neighbors = mapNeighbors(structure_pbc.distances, structure_pbc.symbols, settings.bond_thr, settings.homoatomic_bonding);

    // for (int j = 0; j < 32; ++j)
    //{
    //      int k = 1984+j;
    //      for (int i = 0; i < structure_pbc.neighbors[k].size(); ++i)
    //      {
    //          int n = structure_pbc.neighbors[k][i] ;
    //          std::cout <<
    //          j << " " <<
    //          k << " " <<
    //          structure_pbc.symbols[k] << " " <<
    //          i << " " <<
    //          n << " " <<
    //          structure_pbc.symbols[n] << " " <<
    //          structure_pbc.image[n][0] << " " <<
    //          structure_pbc.image[n][1] << " " <<
    //          structure_pbc.image[n][2] << " " <<
    //          structure_pbc.coordinates[n][0] << " " <<
    //          structure_pbc.coordinates[n][1] << " " <<
    //          structure_pbc.coordinates[n][2] << " " <<
    //          structure_pbc.distances[k][n] << std::endl;
    //      }
    //  }

    return {
        {structure}, {structure_pbc}};
}

// Function to perform the embedding process
Structure embed(Settings settings, Structure structure_pbc)
{
    // Loop over atoms in the PBC structure
    for (int i = 0; i < structure_pbc.nat; ++i)
    {
        // Check if the atom is at the central image
        std::vector<int> image = structure_pbc.image[i];
        if ((image[0] == 0) && (image[1] == 0) && (image[2] == 0))
        {
            // Generate paths for the atom
            std::vector<std::vector<std::vector<int>>> all_paths = generatePaths(structure_pbc.neighbors, i, settings.l_path, settings.l_neis);
            structure_pbc.atoms[i].cycles = all_paths[1];
            structure_pbc.atoms[i].self_returning_paths = all_paths[0];
            structure_pbc.atoms[i].neis = all_paths[2];
        }
    }

    writeEmbeddingData(structure_pbc, settings);
    return structure_pbc;
}

std::pair<std::vector<std::vector<double>>, std::vector<int>> PBC3(std::vector<std::vector<double>> &coordinates,
                                                                   std::vector<std::vector<double>> &lattice_parameters,
                                                                   std::vector<int> &symbols,
                                                                   const int &n)
{
    int nat = coordinates.size();

    // Calculate the expected size of the new_coordinates vector
    int expected_size = nat * (2 * n + 1) * (2 * n + 1) * (2 * n + 1);

    // Preallocate space for the new_coordinates vector
    std::vector<std::vector<double>> new_coordinates(expected_size, std::vector<double>(3));
    std::vector<int> new_symbols(expected_size);

    int index = 0; // Index to directly access elements in new_coordinates

    for (int i = -n; i < n + 1; ++i)
    {
        for (int j = -n; j < n + 1; ++j)
        {
            for (int k = -n; k < n + 1; ++k)
            {
                // std::vector<std::vector<double>> displace_vector = lattice_parameters;
                std::vector<int> image_t = {i, j, k};
                std::vector<double> displace_vector = pbc_d(lattice_parameters, image_t);
                std::vector<std::vector<double>> displaced_coordinates = coordinates;

                for (int row = 0; row < nat; ++row)
                {
                    for (int col = 0; col < 3; ++col)
                    {

                        new_coordinates[index][col] = displaced_coordinates[row][col] + displace_vector[col];
                    }
                    new_symbols[index] = symbols[row];
                    ++index; // Increment the index after each row is processed
                }
            }
        }
    }

    return std::make_pair(new_coordinates, new_symbols);
}

bool compare(double a, double b)
{
    return a > b; // Sort in descending order
}

struct Group
{
    int os;
    double x, y, z;
};

bool compareGroups(const Group &a, const Group &b)
{
    double sumA = std::abs(a.x) + std::abs(a.y) + std::abs(a.z);
    double sumB = std::abs(b.x) + std::abs(b.y) + std::abs(b.z);
    return sumA < sumB;
}

void sortBySmallestXSum(std::vector<double> &vec)
{
    // Convert the vector into groups
    std::vector<Group> groups;
    for (int i = 0; i < vec.size(); i += 4)
    {
        Group g;
        g.os = vec[i];
        g.x = vec[i + 1];
        g.y = vec[i + 2];
        g.z = vec[i + 3];
        groups.push_back(g);
    }

    // Sort the groups based on the sum of absolute values of x, y, and z
    std::sort(groups.begin(), groups.end(), compareGroups);

    // Convert the sorted groups back into the vector
    for (int i = 0; i < groups.size(); ++i)
    {
        vec[i * 4] = groups[i].os;
        vec[i * 4 + 1] = groups[i].x;
        vec[i * 4 + 2] = groups[i].y;
        vec[i * 4 + 3] = groups[i].z;
    }
}

void saveVectorToFile(const std::vector<int> &vec, const std::string &filename)
{
    // Open the file for writing
    std::ofstream file(filename);

    // Check if the file is opened successfully
    if (file.is_open())
    {
        // Write each integer from the vector to the file
        for (int i = 0; i < vec.size(); ++i)
        {
            file << vec[i] << " ";
        }
        file.close();
        std::cout << "Vector saved to file successfully." << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

void saveVectorToBinaryFile(const std::vector<std::vector<double>> &data, const std::string &filename)
{
    std::ofstream file(filename, std::ios::out | std::ios::binary);

    if (!file.is_open())
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    // Write the size of the vector
    size_t rows = data.size();
    file.write(reinterpret_cast<const char *>(&rows), sizeof(size_t));

    // Write the size of each sub-vector and then the data
    for (const auto &row : data)
    {
        size_t cols = row.size();
        file.write(reinterpret_cast<const char *>(&cols), sizeof(size_t));
        file.write(reinterpret_cast<const char *>(row.data()), cols * sizeof(double));
    }

    file.close();
}
// Main function
int main(int argc, char *argv[])
{
    std::string input_name;

    // Command line argument handling
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <input_file_name>, using local file config.in instead\n";
        input_name = "/workspace/shafiro1/post/cycles/bin/config.in";
    }
    else
    {
        input_name = argv[1];
    }

    // Load settings from the input file
    Settings settings = loadSettings(input_name);

    // Embedding process (if enabled in the settings)
    if (settings.embed)
    {
        std::string inputname0 = settings.datapath + settings.input_coordinates;
        std::string inputname1 = settings.datapath + settings.input_lattice_parameters;
        std::string inputname2 = settings.datapath + "cycles_pbc.dat";
        std::string inputname3 = settings.datapath + "self_returning_pbc.dat";
        std::string inputname4 = settings.datapath + "neis.dat";

        std::cout << "Embedding new structure" << std::endl;
        const clock_t embed_time_begin = clock();
        std::vector<Structure> structures = initEmbed(inputname0, inputname1, settings);
        Structure structure = structures[0];
        Structure structure_pbc = structures[1];
        embed(settings, structure_pbc);

        structure_pbc.cycles_full = readPaths(inputname2);
        structure_pbc.self_returning_full = readPaths(inputname3);
        structure_pbc.neis_full = readPaths(inputname4);
        // std::vector<std::vector<int>> neighbors_data = readPaths(inputname4);
        // structure_pbc.cycles_bonds = obtainBonds(structure.cycles);
        // structure_pbc.paths_bonds = obtainBonds(structure.self_returning);
        // writeBondingData(structure_pbc, settings);
        for (int i = 0; i < structure_pbc.neis_full.size(); ++i)
        {
            std::vector<int> cycle = structure_pbc.neis_full[i];
            std::vector<int> cycle_new;
            std::vector<int> only_cyc(cycle.begin() + 3, cycle.end());

            int length = cycle[2];
            int atom0 = cycle[0];
            cycle_new.push_back(atom0);
            cycle_new.push_back(0);
            cycle_new.push_back(length);
            int symbol0 = structure_pbc.symbols[atom0];
            for (int j = 0; j < length; ++j)
            {
                int nei = only_cyc[j];
                int symbols_i = structure_pbc.symbols[nei];

                if (atom0 != nei)
                {
                    cycle_new.push_back(only_cyc[j]);
                    cycle_new[1] += 1;
                }
            }
            structure_pbc.neis.push_back(cycle_new);
        }

        if (settings.type == "cycles")
        {
            for (int i = 0; i < structure_pbc.cycles_full.size(); ++i)
            {
                std::vector<int> cycle = structure_pbc.cycles_full[i];
                int length = cycle[2];
                if (length > 1)
                {
                    std::vector<int> symbols_i(length);
                    std::vector<int> only_cyc(cycle.begin() + 3, cycle.end());
                    for (int j = 0; j < length; ++j)
                    {
                        symbols_i[j] = structure_pbc.symbols[only_cyc[j]];
                    }

                    structure_pbc.cycles.push_back(cycle);
                }
                else
                {
                    structure_pbc.cycles.push_back(cycle);
                }
            }
            writeBinaryFile(settings.datapath + "embed_metadata.bin", structure_pbc.cycles);
        }
        if (settings.type == "paths")
        {
            for (int i = 0; i < structure_pbc.self_returning_full.size(); ++i)
            {
                std::vector<int> cycle = structure_pbc.self_returning_full[i];
                int length = cycle[2];
                if (length > 1)
                {
                    std::vector<int> symbols_i(length);
                    std::vector<int> only_cyc(cycle.begin() + 3, cycle.end());
                    for (int j = 0; j < length; ++j)
                    {
                        symbols_i[j] = structure_pbc.symbols[only_cyc[j]];
                    }
                    structure_pbc.self_returning.push_back(cycle);
                }
                else
                {
                    structure_pbc.self_returning.push_back(cycle);
                }
            }
            writeBinaryFile(settings.datapath + "embed_metadata.bin", structure_pbc.self_returning);
        }
        writeBinaryFile(settings.datapath + "neis_metadata.bin", structure_pbc.neis);

        const clock_t embed_time_end = clock();
        float embed_time = float(embed_time_end - embed_time_begin) / CLOCKS_PER_SEC;
        std::cout << "Embedded (t = " << embed_time << " seconds)" << std::endl;
    }

    // Quick mode (if enabled in the settings)
    if (settings.quick)
    {
        // const clock_t t0 = clock();
        std::vector<std::vector<double>> step_structure = readLatticeParamsCSV(settings.runpath + settings.input_coordinates); // Read the coordinates

        int nat = step_structure.size();
        std::vector<std::vector<double>> relative_coordiates(nat, std::vector<double>(3, 0.0));
        std::vector<int> symbols(nat);
        for (int i = 0; i < nat; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                relative_coordiates[i][j] = step_structure[i][j + 1];
            }
            symbols[i] = step_structure[i][0];
        }

        std::vector<std::vector<double>> lattice_parameters = readLatticeParamsCSV(settings.runpath + settings.input_lattice_parameters); // Read lattice parameters
        std::vector<std::vector<int>> embedding_data = readBinaryFile(settings.datapath + "embed_metadata.bin");
        std::vector<std::vector<int>> neis_data = readBinaryFile(settings.datapath + "neis_metadata.bin");

        // const clock_t t1 = clock();

        std::vector<std::vector<double>> coordinates = matrixMultiply(relative_coordiates, lattice_parameters);
        auto result = PBC3(coordinates, lattice_parameters, symbols, settings.n_pbc);

        std::vector<std::vector<double>> coordinates_pbc = result.first;
        std::vector<int> symbols_pbc = result.second;
        // const clock_t t2 = clock();

        int embedding_size = embedding_data.size();
        int neis_size = neis_data.size();
        int atom0 = neis_data[0][0];

        std::map<std::pair<int, int>, std::vector<double>> emb;
        std::map<std::pair<int, int>, int> emb_count;
        std::vector<int> emb_dimensions(settings.l_path, 0);

        // Loop over embedding data
        for (int i = 0; i < embedding_size; ++i)
        {
            const std::vector<int> &cycle = embedding_data[i];
            const std::vector<int> only_cyc(cycle.begin() + 3, cycle.end());
            const double length = cycle[2];
            const double atom = cycle[0];

            if (length > 1)
            {
                double ds = 0;
                for (int j = 0; j < length - 1; ++j)
                {
                    const std::vector<double> &coords0 = coordinates_pbc[only_cyc[j]];
                    const std::vector<double> &coords1 = coordinates_pbc[only_cyc[j + 1]];
                    const double dx = coords0[0] - coords1[0];
                    const double dy = coords0[1] - coords1[1];
                    const double dz = coords0[2] - coords1[2];
                    ds += sqrt(dx * dx + dy * dy + dz * dz);
                }
                emb[std::make_pair(atom - atom0, length - 1)].push_back(exp(-ds));
            }
        }

        // Sort and populate emb_count
        for (auto &pair : emb)
        {
            std::sort(pair.second.begin(), pair.second.end(), compare);
            const int s = pair.second.size();
            emb_count[pair.first] = s;
            emb_dimensions[pair.first.second] = std::max(emb_dimensions[pair.first.second], s);
        }

        // Populate emb_matrix
        std::vector<std::vector<double>> emb_matrix(nat, std::vector<double>());
        for (int i = 0; i < settings.l_path; ++i)
        {
            std::vector<std::vector<double>> matrix(nat, std::vector<double>(emb_dimensions[i], 0.0));

            // Populate the matrix
            for (const auto &pair : emb)
            {
                if (pair.first.second == i)
                {
                    const int row = pair.first.first;
                    for (size_t j = 0; j < pair.second.size(); ++j)
                    {
                        matrix[row][j] = pair.second[j];
                    }
                }
            }

            // Concatenate horizontally to emb_matrix
            for (int row = 0; row < nat; ++row)
            {
                emb_matrix[row].insert(emb_matrix[row].end(), matrix[row].begin(), matrix[row].end());
            }
        }

        std::vector<int> orders(coordinates_pbc.size(), 1);

        std::map<std::pair<int, int>, std::vector<double>> neis;
        std::map<std::pair<int, int>, std::vector<double>> neis_distance;
        std::map<std::pair<int, int>, std::vector<size_t>> neis_idx;

        for (int i = 0; i < neis_size; ++i)
        {
            std::vector<int> neis_i = neis_data[i];
            std::vector<int> only_neis(neis_i.begin() + 3, neis_i.end());
            double atom = neis_i[0];
            double length = orders[atom];
            double n_neis = neis_i[1];

            orders[atom] += 1;
            std::vector<double> coords_i = coordinates_pbc[atom];

            for (int j = 0; j < n_neis; ++j)
            {
                int idx = only_neis[j];
                std::vector<double> nei_coords = coordinates_pbc[idx];
                std::pair<double, double> im(atom - atom0, length);

                neis[im].push_back(symbols_pbc[idx]);
                neis[im].push_back(nei_coords[0] - coords_i[0]);
                neis[im].push_back(nei_coords[1] - coords_i[1]);
                neis[im].push_back(nei_coords[2] - coords_i[2]);
                if (atom == 1984)
                {
                    std::cout << length << " " << idx << " " << symbols_pbc[idx] << "," << nei_coords[0] << "," << nei_coords[1] << "," << nei_coords[2] << std::endl;
                }
            }
        }
        std::vector<std::vector<double>> nei_matrix(nat, std::vector<double>());
        std::vector<size_t> neis_dimensions(settings.l_neis, 0);

        for (auto &pair : neis)
        {
            neis_dimensions[pair.first.second] = std::max(neis_dimensions[pair.first.second], pair.second.size());
            sortBySmallestXSum(pair.second);
        }

        for (int i = 0; i < settings.l_neis; ++i)
        {
            std::vector<std::vector<double>> matrix(nat, std::vector<double>(neis_dimensions[i], 0.0));
            for (auto &pair : neis)
            {
                if (pair.first.second == i)
                {
                    const int row = pair.first.first;

                    for (size_t j = 0; j < pair.second.size(); ++j)
                    {
                        matrix[row][j] = pair.second[j];
                    }
                }
            }

            // Concatenate horizontally to emb_matrix
            for (int row = 0; row < nat; ++row)
            {
                nei_matrix[row].insert(nei_matrix[row].end(), matrix[row].begin(), matrix[row].end());
            }
        }

        saveVectorToBinaryFile(emb_matrix, settings.runpath + settings.name + "_embed.bin");
        saveVectorToBinaryFile(nei_matrix, settings.runpath + settings.name + "_neis.bin");
    }

    return 0;
}
