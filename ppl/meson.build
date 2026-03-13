gmp_dep = dependency('gmp')
gmpxx_dep = dependency('gmpxx')

inc_cysignals = run_command(
    py,
    [
        '-c',
        '''
from os.path import relpath
import cysignals
path = cysignals.__file__.replace('__init__.py', '')
try:
  print(relpath(path))
except Exception:
  print(path)
    '''.strip(),
    ],
    check: true,
).stdout().strip()
cysignals_dep = declare_dependency(include_directories: inc_cysignals)

inc_gmpy2 = run_command(
    py,
    [
        '-c',
        '''
from os.path import relpath
import gmpy2
path = gmpy2.__file__.replace('__init__.py', '')
try:
  print(relpath(path))
except Exception:
  print(path)
    '''.strip(),
    ],
    check: true,
).stdout().strip()
gmpy2_dep = declare_dependency(include_directories: inc_gmpy2)

ppl_shim_inc = include_directories('.')

common_deps = [ppl_dep, gmp_dep, gmpxx_dep, cysignals_dep, gmpy2_dep]

extensions = {
    'linear_algebra': 'linear_algebra.pyx',
    'mip_problem': 'mip_problem.pyx',
    'polyhedron': 'polyhedron.pyx',
    'generator': 'generator.pyx',
    'constraint': 'constraint.pyx',
    'congruence': 'congruence.pyx',
    'bit_arrays': 'bit_arrays.pyx',
}

foreach ext_name, pyx_file : extensions
    py.extension_module(
        ext_name,
        pyx_file,
        'ppl_shim.cc',
        dependencies: common_deps,
        include_directories: ppl_shim_inc,
        override_options: ['cython_language=cpp'] + ppl_cpp_std_override,
        cython_args: ['-X', 'freethreading_compatible=True'],
        install: true,
        install_dir: py.get_install_dir() / 'ppl',
    )
endforeach

py.install_sources(
    '__init__.py',
    subdir: 'ppl',
)

# Install .pxd files so downstream Cython code can cimport ppl
py.install_sources(
    'ppl_decl.pxd',
    'linear_algebra.pxd',
    'mip_problem.pxd',
    'polyhedron.pxd',
    'generator.pxd',
    'constraint.pxd',
    'congruence.pxd',
    'bit_arrays.pxd',
    'ppl_shim.hh',
    subdir: 'ppl',
)
