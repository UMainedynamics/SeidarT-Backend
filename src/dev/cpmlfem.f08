module cpmlfem 

    use iso_fortran_env, only: real64
    
    implicit none
    
    contains
    
    node_coords(3, nNodes)
    elem_nodes(nLocal, nElems)
    material_id(nElems)
    
    epsilon(nmat, 3, 3) 
    sigma(nmat, 3, 3)
    mu(nmat, 3, 3)
    
    

end module cpmlfem