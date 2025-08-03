from flask_sqlalchemy import SQLAlchemy
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin

db = SQLAlchemy()

class User(UserMixin, db.Model):
    __tablename__ = 'user'
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(150), unique=True, nullable=False)
    password_hash = db.Column(db.String(128), nullable=False)
    email = db.Column(db.String(150), unique=True, nullable=False)

    analyses = db.relationship('GenomeAnalysis', backref='user', lazy=True)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password, method='pbkdf2:sha256')  # método de hash para la contraseña

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

class GenomeAnalysis(db.Model):
    __tablename__ = 'genome_analysis'
    id = db.Column(db.Integer, primary_key=True)
    filename = db.Column(db.String(150))
    length = db.Column(db.Integer)
    count_a = db.Column(db.Integer)
    count_t = db.Column(db.Integer)
    count_c = db.Column(db.Integer)
    count_g = db.Column(db.Integer)
    gc_content = db.Column(db.Float)
    at_gc_ratio = db.Column(db.Float)
    reverse_complement = db.Column(db.Text)
    transcription = db.Column(db.Text)
    translation = db.Column(db.Text)

    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=True)
    timestamp = db.Column(db.DateTime, default=db.func.now())
