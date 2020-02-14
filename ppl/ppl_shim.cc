#include "ppl_shim.hh"

Poly_Gen_Relation* new_relation_with(const Polyhedron &p, const Generator &g)
{
  return new Poly_Gen_Relation(p.relation_with(g));
}

Poly_Con_Relation* new_relation_with(const Polyhedron &p, const Constraint &c)
{
  return new Poly_Con_Relation(p.relation_with(c));
}

Congruence modulo(const Linear_Expression &expr, mpz_class mod)
{
    return (expr %= 0) / mod;
}
