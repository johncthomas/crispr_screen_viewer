from setuptools import setup, find_packages

setup(
    name = 'crispr_screen_viewer',
    version = '0.2.3',
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires = [
        'dash', 'statsmodels','numpy','scipy', 'pandas','pyyaml'
    ],
    python_requires = '>=3.6',
    #scripts=['crispr_tools/crispr_pipeline.py', 'crispr_tools/count_reads.py'],
    include_package_data=False,
)