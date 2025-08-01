//read nodes from grid
std::ifstream gridFileStream(gridFile); // Open the grid file for reading
std::cout<<gridFile<<"\n";
if (!gridFileStream.is_open()) {
    std::cerr << "Unable to open grid file: " << gridFile << std::endl;
    return; // Exit if file cannot be opened
}

// Read nodes from the grid file
std::string line; // To read each line from the CSV
std::getline(gridFileStream, line); // Skip header if present
for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
        Node& node = nodes[i * ny + j];
        std::getline(gridFileStream, line); // Read the entire line
        std::stringstream ss(line); // Create a string stream from the line
        std::string token;

        // Read x and y coordinates, volume, and skewness
        std::getline(ss, token, ','); // Read x_node
        node.x = std::stold(token); // Convert to long double
        std::getline(ss, token, ','); // Read y_node
        node.y = std::stold(token); // Convert to long double
        std::getline(ss, token, ','); // Read volume (not used)
        std::getline(ss, token, ','); // Read skewness (not used)
    }
}
gridFileStream.close(); // Close the grid file after reading