#ifndef PPL_SHIM__H
#define PPL_SHIM__H



#include <ppl.hh>

using namespace Parma_Polyhedra_Library;

// Poly_Gen_Relation/Poly_Con_Relation have no default constructor
Poly_Gen_Relation* new_relation_with(const Polyhedron &p, const Generator &g);
Poly_Con_Relation* new_relation_with(const Polyhedron &p, const Constraint &c);


// Iterator for Generator_System
typedef Generator_System::const_iterator* gs_iterator_ptr;
void delete_gs_iterator(gs_iterator_ptr);

// Iterator for Constraint_System
typedef Constraint_System::const_iterator* cs_iterator_ptr;
void delete_cs_iterator(cs_iterator_ptr);

// Iterator for MIP_Problem
typedef MIP_Problem::const_iterator* mip_cs_iterator_ptr;
void delete_mip_cs_iterator(mip_cs_iterator_ptr mip_csi_ptr);

Constraint_System* mip_constraints(const MIP_Problem &pb);

#endif