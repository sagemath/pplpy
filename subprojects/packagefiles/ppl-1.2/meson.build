project(
    'ppl',
    'cpp',
    version: '1.2',
    meson_version: '>=0.63.0',
)

mod = import('unstable-external_project')

p = mod.add_project(
    'configure',
    configure_options: [
        '--disable-shared',
        '--enable-static',
        '--disable-documentation',
    ],
    env: ['CXXFLAGS=-DNDEBUG -O2 -g0 -fPIC'],
)

ppl_dep = p.dependency('ppl')
meson.override_dependency('ppl', ppl_dep)
