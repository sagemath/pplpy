# ---------- Test environment ----------
# Built ppl extensions live in builddir/ppl/; make them importable
test_env = environment()
test_env.prepend('PYTHONPATH', meson.project_build_root())

# ---------- Pure-Python tests ----------
test(
    'doctests',
    py,
    args: [files('runtests.py')],
    env: test_env,
    workdir: meson.current_source_dir(),
    timeout: 120,
)

test(
    'test_variable',
    py,
    args: ['-m', 'unittest', 'test_variable'],
    env: test_env,
    workdir: meson.current_source_dir(),
)

test(
    'test_constraint',
    py,
    args: ['-m', 'unittest', 'test_constraint'],
    env: test_env,
    workdir: meson.current_source_dir(),
)

# ---------- Cython integration tests ----------
# Include dirs: ../ppl covers both source (for ppl_shim.hh) and build (for
# generated Cython headers).
ppl_src_inc = include_directories('../ppl')

foreach test_name : ['testpplpy', 'testpplpy2']
    ext = py.extension_module(
        test_name,
        test_name + '.pyx',
        dependencies: common_deps,
        include_directories: ppl_src_inc,
        override_options: ['cython_language=cpp'] + ppl_cpp_std_override,
        cython_args: [
            '--directive',
            'language_level=3',
            '-I',
            meson.project_source_root(),
        ],
        install: false,
    )

    test(
        'cython_' + test_name,
        py,
        args: [
            '-c',
            'import sys; sys.path.insert(0, "' + meson.current_build_dir() + '"); ' + 'import ' + test_name + '; ' + test_name + '.test(); ' + test_name + '.example()',
        ],
        env: test_env,
        depends: [ext],
    )
endforeach
