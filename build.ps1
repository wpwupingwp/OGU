rm -Force -Recurse dist
rm -Force -Recurse build
rm -Force -Recurse bf.spec
rm -Force -Recurse *.egg-info
venv/Scripts/Activate.ps1
&pyinstaller OGU/__main__.py -c -F -n bf `
--add-data '.\OGU\data\*.csv;.\OGU\data' `
--add-data '.\OGU\data\*.png;.\OGU\data' `
--add-data '.\venv\Lib\site-packages\primer3\src;.\primer3\src'
&deactivate
&dist/bf.exe
&python -m build -n --wheel
