from setuptools import setup, find_packages

setup(
    name = 'crispr_screen_viewer',
    version = '0.8.0',
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests', 'junk_probably']),
    install_requires = [
        'dash>2', 'numpy','scipy', 'pandas', 'flask', 'plotly',
        'dash-core-components', 'dash-bootstrap-components'
    ],
    scripts = ['crispr_screen_viewer/launch.py'],
    python_requires = '>=3.6',
    include_package_data=False,
)