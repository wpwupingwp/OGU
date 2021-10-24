#!/usr/bin/python3
# database
SQLALCHEMY_DATABASE_URI = 'sqlite:///web.db'
SQLALCHEMY_TRACK_MODIFICATIONS = False
# max filesize 10mb
MAX_CONTENT_LENGTH = 10 * 1024 * 1024
CSRF_ENABLED = True
# safe
SECRET_KEY = '2021'
# bootstrap
BOOTSTRAP_SERVE_LOCAL = True
# database
SQLALCHEMY_DATABASE_URI = 'sqlite:///mai.db'
SQLALCHEMY_TRACK_MODIFICATIONS = False
# upload
UPLOAD_FOLDER = root / 'upload'
UPLOADED_FILE = UPLOAD_FOLDER / 'seq'
for i in UPLOAD_FOLDER, UPLOADED_FILE_DEST:
    if not i.exists():
        i.mkdir()
