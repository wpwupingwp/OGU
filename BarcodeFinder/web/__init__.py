#!/usr/bin/python3

from flask import Flask
from flask_bootstrap import Bootstrap
from flask_login import LoginManager
from pathlib import Path


root = Path.cwd()
app = Flask(__name__)
bootstrap = Bootstrap(app)
lm =LoginManager()
lm.init_app(app)

app.config.from_pyfile('config.py')

from web import config
from web import views
from web.database import db

app.run()