example_env = environment()
example_env.prepend('PYTHONPATH', meson.current_build_dir())

example = py.extension_module(
    'cython_example',
    'cython_example.pyx',
    dependencies: common_deps,
    include_directories: include_directories('../ppl'),
    override_options: ['cython_language=cpp'] + ppl_cpp_std_override,
    cython_args: [
        '--directive',
        'language_level=3',
        '-I',
        meson.project_source_root(),
    ],
    build_by_default: false,
    install: false,
)

test(
    'cython_example',
    py,
    args: ['-m', 'pytest', files('test_cython_example.py')],
    env: example_env,
    depends: [example],
)
