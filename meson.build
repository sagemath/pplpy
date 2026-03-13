project(
    'pplpy',
    'cpp',
    'cython',
    version: '0.9.0',
)

py = import('python').find_installation(pure: false)

# Try to find system PPL via ppl-config script, fall back to bundled subproject
ppl_config = find_program('ppl-config', required: false)
if ppl_config.found()
    ppl_cppflags = run_command(ppl_config, '--cppflags', check: true).stdout().strip().split()
    ppl_prefix = run_command(ppl_config, '--prefix', check: true).stdout().strip()
    ppl_libdir = ppl_prefix / 'lib'
    ppl_dep = declare_dependency(
        compile_args: ppl_cppflags,
        link_args: ['-L' + ppl_libdir, '-lppl'],
    )
    # System PPL 1.2 headers use std::mem_fun_ref, removed in C++17.
    # Compile consumer code with C++11 to stay compatible.
    ppl_cpp_std_override = ['cpp_std=c++11']
else
    ppl_dep = dependency('ppl')
    ppl_cpp_std_override = []
endif

subdir('ppl')
subdir('tests')
