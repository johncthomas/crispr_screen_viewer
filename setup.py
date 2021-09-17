from setuptools import setup, find_packages

setup(
    name = 'crispr_screen_viewer',
    version = '0.4.2',
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests', 'junk_probably']),
    install_requires = [
        'dash==1.21', 'statsmodels','numpy','scipy', 'pandas','pyyaml', 'attrdict'
    ],
    scripts = ['crispr_screen_viewer/launch_msgv.py', 'crispr_screen_viewer/launch_screen_explorer.py'],
    python_requires = '>=3.6',
    include_package_data=False,
)