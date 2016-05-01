from distutils.core import setup, Extension

pymcfiedler = Extension('pymcfiedler',
                        libraries = ['skit', 'csparse', 'mcfiedler'],
                        include_dirs = ['/usr/lib/python2.7/dist-packages/numpy/f2py/src/'],
                        library_dirs = ['../../fortran-pocs/TraceMIN_Fiedler'],
                        sources = ['pymcfiedlermodule.c'])

setup (name = 'MCFiedler',
              version = '1.0',
              description = 'Python interface to TraceMIN-Fiedler',
              ext_modules = [pymcfiedler])
