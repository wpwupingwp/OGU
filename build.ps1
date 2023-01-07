rm -Force -Recurse dist
rm -Force -Recurse build
rm -Force -Recurse bf.spec
rm -Force -Recurse *.egg-info
venv/Scripts/Activate.ps1
&pyinstaller BarcodeFinder/__main__.py -c -F -n bf --add-binary 'E:\Linux\barcodefinder\BarcodeFinder\data;data'
&deactivate
&dist/bf.exe
&python -m build -n --wheel
