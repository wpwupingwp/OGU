$version="3.7","3.8","3.9","3.10","3.11"
$src="R:\build\primer3-py-0.6.1"
# edit MAKE_BIN in setup.py if necessary
# use mingw32-make https://github.com/niXman/mingw-builds-binaries/releases/tag/12.2.0-rt_v10-rev2
$python="python"
&cd ${src}
foreach (${f} in ${version})
{
    &py -${f} -m venv ${f}
    &${f}/Scripts/Activate.ps1
    &pip install -U cython build setuptools wheel
    &${python} -m build -n --wheel
    &deactivate
}
&ls ${src}\dist