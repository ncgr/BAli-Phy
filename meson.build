# https://github.com/mesonbuild/meson/blob/master/docs/markdown/Porting-from-autotools.md
# http://mesonbuild.com/Syntax.html
# http://mesonbuild.com/howtox.html#set-default-cc-language-version
project('bali-phy', 'cpp', version: '3.0-beta5', license: 'GPLv2', default_options : ['cpp_std=c++14'])

cpp = meson.get_compiler('cpp')

add_global_arguments('-DPACKAGE_VERSION="'+meson.project_version()+'"', language : 'cpp')
add_global_arguments('-D_ARCH_="'+host_machine.system()+' ' + host_machine.cpu() + '"', language : 'cpp')
add_global_arguments(['-Wall','-Wextra','-Wno-sign-compare','-Woverloaded-virtual','-Wstrict-aliasing'], language : 'cpp')

if host_machine.system() == 'osx'
  add_global_arguments('-fvisibility=default', language : 'cpp')
endif

eigen = include_directories('include', is_system: true)

cairo = dependency('cairo', required: false)

## we need to build libdl if we are on windows (cygwin or mingw64, and then find the library
libdl = cpp.find_library('dl')

install_subdir('modules', install_dir: 'lib/bali-phy')
install_subdir('functions', install_dir: 'lib/bali-phy')
install_subdir('help', install_dir: 'lib/bali-phy')

# for boost and eigen, we really want something like yes, no, auto* -- default to auto.
# if we have 'auto', then build our own boost if we want stl-debug and bounds_checking
if get_option('with-system-boost')
endif

if get_option('with-system-eigen')
  eigen = include_directories('')
else
  eigen = include_directories('include', is_system: true)
endif

subdir('src')

