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
    submit = m.SubmitField()


class LoginForm(FlaskForm):
    username = m.StringField('Username', validators=[
        v.input_required(), v.email()])
    password = m.PasswordField('Password', validators=[
        v.input_required(), v.length(min=4)])
    submit = m.SubmitField()


class RawCmd(FlaskForm):
    title = m.StringField('Job Title', validators=[v.input_required()])
    module = m.StringField('Module', validators=[v.input_required()])
    command = m.TextAreaField('Command', validators=[v.input_required()])
    submit = m.SubmitField('Submit')


class Gb2fastaForm(FlaskForm):
    genbank_file = m.StringField('GenBank files')
    output = m.StringField('Output path', validators=[v.input_required()])
    no_divide = m.BooleanField('No divide', default=True)
    rename = m.BooleanField('Rename genes', default=False)
    unique = m.BooleanField('Remove repeat sequences', default=False)
    gene = m.StringField('Gene')
    taxon = m.StringField('Taxonomy')
    min_len = m.StringField('Minimum sequence length', default=1)
    max_len = m.StringField('maximum sequence length', default=300000)
    start_date = DateField('Start date')
    end_date = DateField('End date')
    organelle = m.SelectField('Organelle type', choices=(
        'ignore', 'both', 'mitochondrion', 'plastid'))
    refseq = m.SelectField('Use RefSeq database', choices=('both', 'only', 'no'))
    number = m.StringField('Number of records to download', default=0)
    submit = m.SubmitField()


class PostForm(FlaskForm):
    username = m.StringField('Username', validators=[v.input_required('')])
    content = m.TextAreaField('Content', validators=[v.input_required('')])
    secret = m.BooleanField('Secret')
    submit = m.SubmitField()
