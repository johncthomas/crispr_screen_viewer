from setuptools import setup, find_packages

setup(
    name = 'crispr_screen_viewer',
    version = '0.2.4b1',
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires = [
        'dash', 'statsmodels','numpy','scipy', 'pandas','pyyaml', 'attrdict'
    ],
    scripts = ['launch_scripts/launch_volcano.py', 'launch_scripts/launch_jacks_scatter.py'],
    python_requires = '>=3.6',
    #scripts=['crispr_tools/crispr_pipeline.py', 'crispr_tools/count_reads.py'],
    include_package_data=False,
)