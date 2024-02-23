# currently powershell should use 7.2
cd E:\Linux\barcodefinder
python -m venv build_win
./build_win/scripts/activate.ps1
#pip install -r requirements.txt
#pip install pyinstaller altgraph
cd src
&pyinstaller -n OGU --add-data 'OGU/data/*.csv:OGU/data' --clean -F --hidden-import OGU OGU/__main__.py
deactivate
cd ..
./dist/OGU.exe