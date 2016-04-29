from numpy.distutils.core import Extension, setup

pymcfiedler = Extension('pymcfiedler',
                        libraries = ['mcfiedler'],
                        include_dirs = ['/usr/lib/python2.7/dist-packages/numpy/f2py/src/'],
                        library_dirs = [ ],
                        sources = ['pymcfiedlermodule.c'])

setup (name = 'MCFiedler',
              version = '1.0',
              description = 'Python interface to TraceMIN-Fiedler',
              ext_modules = [pymcfiedler])
