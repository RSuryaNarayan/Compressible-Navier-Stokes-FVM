// loop over faces and push the face into the cells to the left and right
for (int iif=0; iif<nfaces; iif++)
{
    Face face = faces[iif];
    int left = face.i_L;
    int right = face.i_R;
    if (left>=0) cells[left].faces_index.push_back(face.index);
    if (right>=0) cells[right].faces_index.push_back(face.index); 

    //check bc assignment
    // int bc_id = face.bc_index;
    // if (bc_id==4) std::cout<<"Face "<<face.index<<"\n";
}