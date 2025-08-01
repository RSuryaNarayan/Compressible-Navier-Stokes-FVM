for (int i = 0; i < nx - 1; ++i) {
    for (int j = 0; j < ny - 1; ++j) {

        // Populate Cell Information 
        // cell centers
        Cell cell;
        cell.x_c = (nodes[i * ny + j].x + nodes[i * ny + (j + 1)].x + 
                        nodes[(i + 1) * ny + j].x + nodes[(i + 1) * ny + (j + 1)].x) / 4.0;
        cell.y_c = (nodes[i * ny + j].y + nodes[i * ny + (j + 1)].y + 
                        nodes[(i + 1) * ny + j].y + nodes[(i + 1) * ny + (j + 1)].y) / 4.0;
        
        // Store nodes in clockwise order
        cell.nodes = {
            &nodes[i * ny + j],         // Bottom-left
            &nodes[i * ny + (j + 1)],   // Bottom-right
            &nodes[(i + 1) * ny + (j + 1)], // Top-right
            &nodes[(i + 1) * ny + j]    // Top-left
        };

        // Calculate the area of the cell using the shoelace formula
        long double area = 0.5 * std::abs(
            nodes[i * ny + j].x * nodes[i * ny + (j + 1)].y +
            nodes[i * ny + (j + 1)].x * nodes[(i + 1) * ny + (j + 1)].y +
            nodes[(i + 1) * ny + (j + 1)].x * nodes[(i + 1) * ny + j].y +
            nodes[(i + 1) * ny + j].x * nodes[i * ny + j].y -
            (nodes[i * ny + (j + 1)].x * nodes[i * ny + j].y +
                nodes[(i + 1) * ny + (j + 1)].x * nodes[i * ny + (j + 1)].y +
                nodes[(i + 1) * ny + j].x * nodes[(i + 1) * ny + (j + 1)].y +
                nodes[i * ny + j].x * nodes[(i + 1) * ny + j].y)
        );

        cell.volume = area; // Assign the calculated area to the cell's volume

        // Calculate unit normals for each face
        long double sum_nxS = 0;
        long double sum_nyS = 0;
        std::vector<long double> normals; // Store normals for the cell faces
        for (int k = 0; k < 4; ++k) { // 4 faces
            long double x1 = cell.nodes[k]->x;
            long double y1 = cell.nodes[k]->y;
            long double x2 = cell.nodes[(k + 1) % 4]->x; // Next node (wrap around)
            long double y2 = cell.nodes[(k + 1) % 4]->y;

            long double nx = -(y2 - y1); // Normal x-component
            long double ny = (x2 - x1);  // Normal y-component
            long double magnitude = sqrt(nx * nx + ny * ny); // Magnitude
            
            cell.S.push_back(magnitude); // Push face area into cell.S

            // Normalize the normal vector
            if (magnitude != 0) {
                nx /= magnitude;
                ny /= magnitude;
            }

            cell.n_x.push_back(nx);
            cell.n_y.push_back(ny);

            sum_nxS = sum_nxS + nx*magnitude;
            sum_nyS = sum_nyS + ny*magnitude;
        }
        
        cell.faceVectorSum = (sqrt(sum_nxS*sum_nxS + sum_nyS*sum_nyS)); //pushback faceVectorSum

        // Populate index of the cell and its neighbors
        cell.index = i * (ny - 1) + j; // Set the index of the current cell

        // Populate neighbors
        // if (i > 0) cell.neighbors.push_back((i - 1) * (ny - 1) + j); // Top neighbor
        // if (i < nx - 2) cell.neighbors.push_back((i + 1) * (ny - 1) + j); // Bottom neighbor
        // if (j > 0) cell.neighbors.push_back(i * (ny - 1) + (j - 1)); // Left neighbor
        // if (j < ny - 2) cell.neighbors.push_back(i * (ny - 1) + (j + 1)); // Right neighbor

        // Print unit normals for each face of the cell
        // std::cout << "Unit Normals: ";
        // for (int face = 0; face < 4; ++face) {
        //     std::cout << "(" << cell.n_x[face] << ", " << cell.n_y[face] << ") ";
        // }
        // std::cout << std::endl; // New line after printing normals

        // // Print neighbors for each cell
        // std::cout << "Cell index: " << cell.index 
        //           << " Location: (" << cell.x_c << ", " << cell.y_c << ") "
        //           << "Neighbors: ";
        // for (const auto& neighbor : cell.neighbors) {
        //     std::cout << neighbor << " ";
        // }
        // std::cout << std::endl; // New line after printing neighbors

        cells.push_back(cell);
    }
}
