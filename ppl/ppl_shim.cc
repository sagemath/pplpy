#include "ppl_shim.hh"

Poly_Gen_Relation* new_relation_with(const Polyhedron &p, const Generator &g)
{
  return new Poly_Gen_Relation(p.relation_with(g));
}

Poly_Con_Relation* new_relation_with(const Polyhedron &p, const Constraint &c)
{
  return new Poly_Con_Relation(p.relation_with(c));
}

/************************************************************/

Constraint_System* mip_constraints(const MIP_Problem &pb)
{
	Constraint_System *cs = new Constraint_System();
	for(MIP_Problem::const_iterator it = pb.constraints_begin();
		it != pb.constraints_end();
		++it)
		cs->insert(*it);
	return cs;
}
