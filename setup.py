import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='csrstat',
    version='0.1',
    description='Package for introductory analysis of Spatial Point Processes,\
                 including basic utility functions, summary statistics and\
                 metrics for 2D point processes. To be primarily used to\
                 test the hypothesis of CSR.',
    url='https://github.com/chancehaycock/csrstat',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Chance Haycock, Edward Thorpe-Woods',
    author_email='chancehaycock@me.com, t_thorpewoods@hotmail.co.uk',
    license='MIT',
    packages=setuptools.find_packages(),
    classifiers=["Programming Language :: Python :: 3",
                 "License :: OSI Approved :: MIT License",
                 "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
      'pointpats==2.1.0',
      'numpy==1.18.4',                                                                   
      'scipy==1.10.0',                                                                    
      'matplotlib==3.2.1',                                                               
      'pytest==5.4.3',                                                                   
      'pandas==1.0.4',                                                                
      'seaborn==0.10.1',
      'utm==0.5.0'
    ] 
)
