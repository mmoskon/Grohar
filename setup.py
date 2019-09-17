#https://www.smallsurething.com/a-really-simple-guide-to-packaging-your-pyqt-application-with-cx_freeze/

from cx_Freeze import setup, Executable
 
import os

os.environ['TCL_LIBRARY'] = r"C:\ProgramData\Anaconda3\tcl\tcl8.6"
os.environ['TK_LIBRARY'] = r"C:\ProgramData\Anaconda3\tcl\tk8.6"

# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = dict(packages = [], excludes = [])
 
import sys
base = 'Win32GUI' if sys.platform=='win32' else None
 
executables = [
    Executable('top.py', base=base)
]
 
setup(
    name='top',
    version = '0.1',
    description = 'Metabolic Vizualisation Tool',
    options = dict(build_exe = buildOptions),
    executables = executables
)