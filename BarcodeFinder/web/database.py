#!/usr/bin/python3

from datetime import datetime
from flask_sqlalchemy import SQLAlchemy

from web import app
from flask_wtf import FlaskForm

import flask_login as fl

db = SQLAlchemy(app)


class User(db.Model, fl.UserMixin):
    __tablename__ = 'user'
    user_id = db.Column(db.Integer, primary_key=True)
    # email
    username = db.Column(db.String(100), unique=True)
    email = db.Column(db.String(100), unique=True)
    password = db.Column(db.String(100))
    register_date = db.Column(db.DateTime)
    address = db.Column(db.String(100))
    jober_id = db.relationship('Command', backref='user')

    def __init__(self, username, password, address=''):
        self.username = username
        self.password = password
        self.register_date = datetime.utcnow()
        self.address = address

    def __repr__(self):
        return f'{self.username}'

    def get_id(self):
        return str(self.user_id)


class Command(db.Model):
    __tablename__ = 'command'
    id = db.Column(db.Integer, primary_key=True)
    title = db.Column(db.String(100), default='Unknown', nullable=False)
    user_id = db.Column(db.ForeignKey('user.user_id'))
    # 1, 12, 123, 23, 2, 3
    module = db.Column(db.Integer, default=123, nullable=False)
    command = db.Column(db.Text, nullable=False)
    date = db.Column(db.DateTime)

    def __str__(self):
        return f'{self.command_id}: {self.title}'

    @staticmethod
    def from_form(form, user_id):
        cmd = Command()
        if isinstance(form, FlaskForm):
            cmd.user_id = user_id
            cmd.title = 'test'
            cmd.module = 1
            cmd.command = str(dict(form.data))
            cmd.date = datetime.now()
        else:
            pass
        return cmd


class Post(db.Model):
    __tablename__ = 'post'
    post_id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(100), nullable=False)
    secret = db.Column(db.Boolean, default=False)
    content = db.Column(db.Text, nullable=False)
    date = db.Column(db.DateTime)

    def __init__(self, username, content, secret):
        self.username = username
        self.content = content
        self.secret = secret
        self.date = datetime.now()

    def __repr__(self):
        return f'{self.post_id}'