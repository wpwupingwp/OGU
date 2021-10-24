#!/usr/bin/python3

from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
import wtforms as m
from wtforms import validators as v
from wtforms.fields.html5 import DateField


class UserForm(FlaskForm):
    username = m.StringField('Username', validators=[v.input_required()])
    email = m.StringField('Email', validators=[
        v.input_required(), v.email()])
    password = m.PasswordField('Password', validators=[
        v.input_required(), v.length(min=4)])
    password2 = m.PasswordField('Password Again', validators=[
        v.input_required(), v.equal_to('password'), v.length(min=4)])
    address = m.StringField('Address', validators=[v.input_required()])
    submit = m.SubmitField('Submit')


class LoginForm(FlaskForm):
    username = m.StringField('Username', validators=[
        v.input_required(), v.email()])
    password = m.PasswordField('Password', validators=[
        v.input_required(), v.length(min=4)])
    submit = m.SubmitField('Submit')


class RawCmd(FlaskForm):
    title = m.StringField('Job Title', validators=[v.input_requried()])
    module = m.StringField('Module', validators=[v.input_requried()])
    command = m.TextAreaField('Command', validators=[v.input_required()])
    submit = m.SubmitField('Submit')


class PostForm(FlaskForm):
    username = m.StringField('Username', validators=[v.input_required('')])
    content = m.TextAreaField('Content', validators=[v.input_required('')])
    secret = m.BooleanField('Secret')
    submit = m.SubmitField('Submit')
