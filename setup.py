import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
        if cmake_version < LooseVersion('3.5.0'):
            raise RuntimeError("CMake >= 3.5.0 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        build_type = os.environ.get("BUILD_TYPE", "Release")
        build_args = ['--config', build_type]

        # Pile all .so in one place and use $ORIGIN as RPATH
 
        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(build_type.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + build_type]
            build_args += ['--', '-j4']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake',
                               '--build', '.',
                               '--target', os.path.basename(ext.name)
                               ] + build_args,
                              cwd=self.build_temp)


setup(
        name='pyMBA',
        version='0.0.1',
        description='Multilevel B-spline interpolation',
        author='Lionel Untereiner',
        author_email='lionel.untereiner@geosiris.com',
        license='MIT',
        url='https://github.com/untereiner/MBA',
        include_package_data=True,
        ext_modules=[CMakeExtension('pyMBA', 'src')],
        cmdclass=dict(build_ext=CMakeBuild),
        zip_safe=False,
        # [
        #     Extension('MBA', ['python/pymba.cpp'],
        #         include_dirs=['python/pybind11/include/', 'include/'],
        #         extra_compile_args=['-O3', '-std=c++11']
        #         )
        #     ]
)

# from setuptools import setup
# import subprocess
# import os
# import sys
# import platform
# from pybind11.setup_helpers import Pybind11Extension, build_ext
# from pybind11 import get_cmake_dir



# class cmake_build_ext(build_ext):
#     def build_extensions(self):
#         # Ensure that CMake is present and working
#         try:
#             out = subprocess.check_output(['cmake', '--version'])
#         except OSError:
#             raise RuntimeError('Cannot find CMake executable')

#         for ext in self.extensions:

#             extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
#             cfg = 'Debug' #if options['--debug'] == 'ON' else 'Release'

#             cmake_args = [
#                 '-DCMAKE_BUILD_TYPE=%s' % cfg,
#                 # Ask CMake to place the resulting library in the directory
#                 # containing the extension
#                 '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
#                 # Other intermediate static libraries are placed in a
#                 # temporary build directory instead
#                 '-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), self.build_temp),
#                 # Hint CMake to use the same Python executable that
#                 # is launching the build, prevents possible mismatching if
#                 # multiple versions of Python are installed
#                 '-DPYTHON_EXECUTABLE={}'.format(sys.executable),
#                 # Add other project-specific CMake arguments if needed
#                 # ...
#                 '-Dmba_ROOT=/home/untereiner/Documents/Developments/geosiris/gtm/v2/MBA/build/cmake/mba'
#             ]

#             # We can handle some platform-specific settings at our discretion
#             if platform.system() == 'Windows':
#                 plat = ('x64' if platform.architecture()[0] == '64bit' else 'Win32')
#                 cmake_args += [
#                     # These options are likely to be needed under Windows
#                     '-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE',
#                     '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
#                 ]
#                 # Assuming that Visual Studio and MinGW are supported compilers
#                 if self.compiler.compiler_type == 'msvc':
#                     cmake_args += [
#                         '-DCMAKE_GENERATOR_PLATFORM=%s' % plat,
#                     ]
#                 else:
#                     cmake_args += [
#                         '-G', 'MinGW Makefiles',
#                     ]

#             # cmake_args += cmake_cmd_args

#             if not os.path.exists(self.build_temp):
#                 os.makedirs(self.build_temp)

#             # Config
#             subprocess.check_call(['cmake'] + cmake_args,
#                                   cwd=self.build_temp)

#             # Build
#             subprocess.check_call(['cmake', '--build', '.', '--config', cfg],
#                                   cwd=self.build_temp)


# __version__ = "0.0.1"

# ext_modules = [
#     Pybind11Extension("pyMBA",
#         ["src/pyMBA.cpp"], include_dirs=["include"],
#         # Example: passing in the version to the compiled code
#         define_macros = [('VERSION_INFO', __version__)],
#         ),
# ]

# setup(
#         name='MBA',
#         version=__version__,
#         description='Multilevel B-spline interpolation',
#         author='Lionel Untereiner',
#         author_email='lionel.untereiner@geosiris.com',
#         license='MIT',
#         url='https://github.com/untereiner/MBA',
#         include_package_data=True,
#         ext_modules=ext_modules,
#         cmdclass={"build_ext": cmake_build_ext},
#         # cmdclass=dict(build_ext=ExtensionBuilder),
#         zip_safe=False
# )