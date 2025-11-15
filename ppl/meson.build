py.install_sources(
  '__init__.py',
  'bit_arrays.pxd',
  'congruence.pxd',
  'constraint.pxd',
  'generator.pxd',
  'linear_algebra.pxd',
  'mip_problem.pxd',
  'polyhedron.pxd',
  'ppl_decl.pxd',
  'ppl_shim.hh',
  subdir: 'ppl'
)

# Shared C++ source file
ppl_shim_lib = static_library(
  'ppl_shim',
  'ppl_shim.cc',
  dependencies: [ppl, gmp],
  cpp_args: [],
)

extension_data = {
  'linear_algebra': files('linear_algebra.pyx'),
  'mip_problem': files('mip_problem.pyx'),
  'polyhedron': files('polyhedron.pyx'),
  'generator': files('generator.pyx'),
  'constraint': files('constraint.pyx'),
  'congruence': files('congruence.pyx'),
  'bit_arrays': files('bit_arrays.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'ppl',
    install: true,
    dependencies: [cysignals, gmpy2_dep, ppl, gmp],
    link_with: ppl_shim_lib,
    override_options : ['cython_language=cpp'],
  )
endforeach
