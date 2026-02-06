project('pplpy',
  ['c', 'cpp', 'cython'],
  version: files('VERSION'),
  license: 'GPL v3',
  default_options: ['cpp_std=c++11', 'python.install_env=auto'],
  meson_version: '>=1.2',
)

# Python module
# https://mesonbuild.com/Python-module.html
py = import('python').find_installation(pure: false)

# Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cython = meson.get_compiler('cython')
# Workaround as described in https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#arithmetic-methods
add_project_arguments('-X c_api_binop_methods=True', language: 'cython')

# Dependencies
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
cysignals = declare_dependency(include_directories: inc_cysignals)

# Find gmpy2 include directory
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

# Find PPL library
ppl = cpp.find_library('ppl', required: true)

# Find GMP library
gmp = cc.find_library('gmp', required: true)

subdir('ppl')
