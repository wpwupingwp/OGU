from web.database import db, User, Post

try:
    db.create_all()
except:
    pass

if User.query.filter_by(username='admin').first() is None:
    u = User('admin', '123456', '')
    db.session.add(u)
if User.query.filter_by(username='guest').first() is None:
    u = User('guest', '123456', '')
    db.session.add(u)
for i in range(5):
    p = Post('admin', f'测试-{i}', False)
    db.session.add(p)
db.session.commit()