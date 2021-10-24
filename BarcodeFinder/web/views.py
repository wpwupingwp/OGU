#!/usr/bin/python3

import flask as f

from class_demo import app
from class_demo.database import Post, db
from class_demo.form import PostForm


@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return f.send_from_directory(app.config['UPLOADED_FILE_DEST'], filename)


@app.route('/')
@app.route('/index')
def index():
    return f.render_template('index.html')


@app.route('/gb2fasta')
def gb2fasta():
    form = Gb2fastaForm()
    if form.validate_on_submit():
        cmd = Command()
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('form.html', form=form)


@app.route('/evaluate')
def gb2fasta():
    form = EvaluateForm()
    if form.validate_on_submit():
        cmd = Command()
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('form.html', form=form)


@app.route('/primer')
def gb2fasta():
    form = PrimerForm()
    if form.validate_on_submit():
        cmd = Command()
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('form.html', form=form)


@app.route('/combine')
def gb2fasta():
    form = CommandForm()
    if form.validate_on_submit():
        cmd = Command()
        db.session.add(cmd)
        db.session.commit()
        f.flash('Submit OK')
        return f.redirect('/')
    return f.render_template('form.html', form=form)


@app.route('/job/<int:page>')
def my_goods(user_id, page=1):
    per_page = 5
    pagination = Command.query.with_entities(
        Command.id, Command.title, Command.user_id, Command.date
    ).order_by(Command.date.desc()).paginate(page=page, per_page=per_page)
    return f.render_template('job.html', pagination=pagination)


@app.route('/post', methods=('POST', 'GET'))
@app.route('/post/<int:page>', methods=('POST', 'GET'))
def index(page=1):
    postform = PostForm()
    pagination = Post.query.order_by(Post.date.desc()).paginate(page=page,
                                                                per_page=5)
    if postform.validate_on_submit():
        post = Post(postform.username.data, postform.content.data,
                    postform.secret.data)
        db.session.add(post)
        db.session.commit()
        f.flash('留言成功')
        return f.redirect('/')
    return f.render_template('index.html', form=postform, pagination=pagination)