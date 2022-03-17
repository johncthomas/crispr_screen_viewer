from setuptools import setup, find_packages

setup(
    name = 'crispr_screen_viewer',
    version = '0.6.0',
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests', 'junk_probably']),
    #todo review requires, do we need attrdict/scipy
    install_requires = [
        'dash>2', 'statsmodels','numpy','scipy', 'pandas', 'flask', 'plotly',
        'dash-html-components', 'dash-core-components'
    ],
    scripts = ['crispr_screen_viewer/main_page.py'],
    python_requires = '>=3.6',
    include_package_data=False,
)