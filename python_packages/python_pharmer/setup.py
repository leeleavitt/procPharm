import setuptools 

setuptools.setup(name='python_pharmer',
      version='1',
      description='Extracting ND2 video data',
      author='Lee Leavitt, Rishi Alluri',
      author_email='lee.leavitt.u@gmail.com',
      packages=setuptools.find_packages(),
      package_data={'':['models/*.h5']},
      install_requires=[
            'pims==0.4.1', 
            'tensorflow==2.3.1', 
            'keras==2.3.1', 
            'pandas==1.0.1', 
            'scikit-image==0.16.2', 
            'numpy==1.18.1',
            'xmltodict==0.12.0'
      ]
)

