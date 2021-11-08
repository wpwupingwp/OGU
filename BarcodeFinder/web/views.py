#!/usr/bin/python3

import flask as f
import flask_login as fl

from web import app, lm
from web.database import Command, Post, db
from web.form import *


@lm.user_loader
def load_user(user_id):
    user = User.query.get(user_id)
    return user


@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return f.send_from_directory(app.config['UPLOADED_FILE_DEST'], filename)


@app.route('/')
@app.route('/index')
def index():
    return f.render_template('index.html')


@app.route('/gb2fasta', methods=('POST', 'GET'))
def gb2fasta():
    form = Gb2fastaForm()
    if form.validate_on_submit():
        cmd = Command.from_form(form, 1)
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('base.html', form=form)


@app.route('/evaluate', methods=('POST', 'GET'))
def evaluate():
    form = EvaluateForm()
    if form.validate_on_submit():
        cmd = Command()
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('base.html', form=form)


@app.route('/primer', methods=('POST', 'GET'))
def primer():
    form = PrimerForm()
    if form.validate_on_submit():
        cmd = Command.from_form(form, 4)
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('base.html', form=form)


@app.route('/combine', methods=('POST', 'GET'))
def combine():
    form = RawCmd()
    if form.validate_on_submit():
        cmd = Command()
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('base.html', form=form)


@app.route('/job/<int:page>')
@app.route('/job/')
def my_goods(page=1):
    per_page = 5
    pagination = Command.query.with_entities(
        Command.id, Command.title, Command.module, Command.user_id, Command.date
    ).order_by(Command.date.desc()).paginate(page=page, per_page=per_page)
    return f.render_template('job.html', pagination=pagination)


@app.route('/post', methods=('POST', 'GET'))
@app.route('/post/<int:page>', methods=('POST', 'GET'))
def post(page=1):
    postform = PostForm()
    pagination = Post.query.order_by(Post.date.desc()).paginate(page=page,
                                                                per_page=5)
    if postform.validate_on_submit():
        post = Post(postform.username.data, postform.content.data,
                    postform.secret.data)
        db.session.add(post)
        db.session.commit()
        f.flash('Post OK.')
        return f.redirect('/post')
    return f.render_template('post.html', form=postform, pagination=pagination)