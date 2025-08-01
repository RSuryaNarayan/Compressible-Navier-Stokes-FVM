// First populate the "y-"faces
int count = 0; 
for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny-1; ++j) {
        Face face;
        face.index = count;
        face.nodes = {
            &nodes[i * ny + j + 1 ], 
            &nodes[i * ny + j + 0 ]
        };
        long double x1 = face.nodes[0]->x;
        long double y1 = face.nodes[0]->y;
        long double x2 = face.nodes[1]->x; 
        long double y2 = face.nodes[1]->y;

        long double n_x = -(y2 - y1); // Normal x-component
        long double n_y = (x2 - x1);  // Normal y-component
        long double magnitude = std::sqrt(n_x * n_x + n_y * n_y); // Magnitude

        // Normalize the normal vector
        if (magnitude != 0) {
            n_x /= magnitude;
            n_y /= magnitude;
        }

        face.n_x = n_x;
        face.n_y = n_y;
        face.S_f = magnitude;

        // Set boundary conditions for the edges of the domain
        if (i == 0) face.bc_index = 1; // Left boundary
        else if (i == nx - 1) face.bc_index = 2; // Right boundary
        else face.bc_index = 0; // Internal face

                        // set left and right indices
        if (face.bc_index==0) // interior cell
        {
            face.i_L = face.index - (ny-1);
            face.i_R = face.index;
        }
        else if (face.bc_index==1) // inlet
        {
            face.i_L = -1; //junk
            face.i_R = face.index;
        }
        else if (face.bc_index==2) // outlet
        {
            face.i_L = face.index - (ny-1);
            face.i_R = -1; // junk
        }

        // std::cout<<"Face "<<face.index<<" cell left: "<<face.i_L<<" x_c_L="<<cells[face.i_L].x_c<<"  cell right: "<<face.i_R<<"\n";
        // Push to the global face vector
        faces.push_back(face);
        count++;
    }
}

// Next populate "x-"faces
for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx-1; ++i) {
        Face face;
        face.index = count;
        face.nodes = {
            &nodes[ (i+0) * ny + j ], 
            &nodes[ (i+1) * ny + j ]
        };
        long double x1 = face.nodes[0]->x;
        long double y1 = face.nodes[0]->y;
        long double x2 = face.nodes[1]->x; 
        long double y2 = face.nodes[1]->y;

        long double n_x = -(y2 - y1); // Normal x-component
        long double n_y = (x2 - x1);  // Normal y-component
        long double magnitude = sqrt(n_x * n_x + n_y * n_y); // Magnitude

        // Normalize the normal vector
        if (magnitude != 0) {
            n_x /= magnitude;
            n_y /= magnitude;
        }

        face.n_x = n_x;
        face.n_y = n_y;
        face.S_f = magnitude;

        // Set boundary conditions for the edges of the domain
        if (j == 0) face.bc_index = 3; // Bottom boundary
        else if (j == ny - 1) face.bc_index = 4; // Top boundary
        else face.bc_index = 0; // Internal face

        // set left and right for face
        if (face.bc_index==0) // interior
        {
            face.i_L = i*(ny-1)+j-1;
            face.i_R = i*(ny-1)+j;
        }
        else if (face.bc_index==3) // bottom boundary
        {
            face.i_L = -1;
            face.i_R = i*(ny-1)+j;
        }
        else if (face.bc_index==4) // top boundary
        {
            face.i_L = i*(ny-1)+j-1;
            face.i_R = -1;
        }

        // Push to the global face vector
        faces.push_back(face);

        //check by printing
        // std::cout<<"Face "<<faces[count].index<<" node 1: ("<<faces[count].nodes[0]->x<<","<<faces[count].nodes[0]->y
        //                                     <<"), node 2: ("<<faces[count].nodes[1]->x<<","<<faces[count].nodes[1]->y<<")\n";
        // std::cout<<"Face "<<face.index<<" cell left: "<<face.i_L<<" x_c_L="<<cells[face.i_L].x_c<<"  cell right: "<<face.i_R<<"\n";

        count++;
    }
}

//update number of faces
nfaces = faces.size();
std::cout<<"\n Number of faces = "<<nfaces;
if (nfaces==2*nx*ny-nx-ny) std::cout<<"\n Number of faces successfully identified!\n";
else {
    std::cout<<"\n something wrong with faces! aborting!!\n";                     
    std::abort(); // Abort the program
}

//update number of cells
ncells = cells.size();
std::cout<<"\n Number of cells = "<<ncells;
if (ncells==(nx-1)*(ny-1)) std::cout<<"\nNumber of cells successfully identified\n";
else {
    std::cout<<"\n something wrong with cells! aborting!!\n";                     
    std::abort(); // Abort the program
}
