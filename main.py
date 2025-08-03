import os 
from flask import Flask, render_template, request, redirect, url_for, send_file, flash
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager, login_user, logout_user, login_required, current_user
from models import db, User, GenomeAnalysis 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from werkzeug.security import check_password_hash, generate_password_hash


UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)
os.makedirs("database", exist_ok=True)


app = Flask(__name__) 
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = os.environ.get('SECRET_KEY', 'dev_secret_key_for_testing_only')


# DATABASE CONFIGURATION
basedir = os.path.abspath(os.path.dirname(__file__))
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///D:/genomeapp/database/genome_analysis.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Initialize the database with the app
db.init_app(app)

# LOGIN CONFIGURATION
login_manager = LoginManager()
login_manager.login_view = 'home' 
login_manager.init_app(app)

@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

def analyze_fasta(file_path):
    record = SeqIO.read(file_path, "fasta")
    seq = record.seq.upper()

    a_count = seq.count("A")
    c_count = seq.count("C")
    t_count = seq.count("T")
    g_count = seq.count("G")
    length = len(seq)
    at = a_count + t_count
    gc = g_count + c_count
    at_gc_ratio = at / gc if gc != 0 else float("inf")
    gc_content = (gc / length) * 100 if length > 0 else 0

    rev_comp = seq.reverse_complement()
    rna = seq.transcribe()
    
    # Define paths
    rev_path = os.path.join("results", "reverse_complement.fasta")
    rna_path = os.path.join("results", "transcription.fasta")

    # Create SeqRecord objects with headers
    rev_record = SeqRecord(Seq(str(rev_comp)), id="reverse_complement", description="")
    rna_record = SeqRecord(Seq(str(rna)), id="transcription", description="")

    # Write valid FASTA files
    SeqIO.write(rev_record, rev_path, "fasta")
    SeqIO.write(rna_record, rna_path, "fasta")

    return {
        "A": a_count,
        "C": c_count,
        "T": t_count,
        "G": g_count,
        "Length": length,
        "AT/GC Ratio": round(at_gc_ratio, 2),
        "GC%": round(gc_content, 2),
        "Reverse complement path": rev_path,
        "Transcription path": rna_path,
        "Translation": str(seq.translate(to_stop=True))
    }
    

@app.route('/', methods=['GET', 'POST'])
def home():
    # Default data cards (used for both GET and no file uploaded)
    default_data_cards = [
        {"title": "A", "value": "Count of <br>Adenine bases", "color": "#DBE9EE"},
        {"title": "C", "value": "Count of <br>Cytosine bases", "color": "#DBE9EE"},
        {"title": "T", "value": "Count of <br>Thymine bases", "color": "#DBE9EE"},
        {"title": "G", "value": "Count of <br>Guanine bases", "color": "#DBE9EE"},
        {"title": "Length", "value": "Total sequence length", "color": "#C0D6DF"},
        {"title": "AT/GC Ratio", "value": "Ratio of AT to GC content", "color": "#C0D6DF"},
        {"title": "GC%", "value": "Percentage of GC content", "color": "#C0D6DF"},
        {"title": "Reverse complement", "value": "Reverse-complement <br>of the sequence", "color": "#4F6D7A"},
        {"title": "Transcription", "value": "RNA sequence from transcription", "color": "#4F6D7A"},
        {"title": "Translation", "value": "Protein sequence", "color": "#4F6D7A"},
    ]
    
    if request.method == 'POST':
        uploaded_file = request.files['genome_file']
        if uploaded_file.filename != '':
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], uploaded_file.filename)
            uploaded_file.save(file_path)

            result = analyze_fasta(file_path)
            analysis = GenomeAnalysis(
                filename=uploaded_file.filename,
                length=result["Length"],
                count_a=result["A"],
                count_t=result["T"],
                count_c=result["C"],
                count_g=result["G"],
                gc_content=result["GC%"],
                at_gc_ratio=result["AT/GC Ratio"],
                reverse_complement=open(result["Reverse complement path"]).read(),
                transcription=open(result["Transcription path"]).read(),
                translation=result["Translation"],
                user_id=None  # Si no has implementado login todavía
                )

            db.session.add(analysis)
            db.session.commit()

            data_cards = [
                {"title": "A", "value": result["A"], "unit": "bases", "color": "#DBE9EE"},
                {"title": "C", "value": result["C"], "unit": "bases", "color": "#DBE9EE"},
                {"title": "T", "value": result["T"], "unit": "bases", "color": "#DBE9EE"},
                {"title": "G", "value": result["G"], "unit": "bases", "color": "#DBE9EE"},
                {"title": "Length", "value": result["Length"], "unit": "bases", "color": "#C0D6DF"},
                {"title": "AT/GC Ratio", "value": result["AT/GC Ratio"], "unit": "ratio", "color": "#C0D6DF"},
                {"title": "GC%", "value": result["GC%"], "unit": "%", "color": "#C0D6DF"},
                {"title": "Reverse complement", "value": "Download", "color": "#4F6D7A", "url": 'reverse_complement.fasta',"is_button": True},
                {"title": "Transcription", "value": "Download", "color": "#4F6D7A", "url": 'transcription.fasta',"is_button": True},
                {"title": "Translation", "value": result["Translation"], "color": "#4F6D7A"},
            ]
            return render_template("index.html", data_cards=data_cards)
        
    # If it's a GET or POST with no file — return default descriptions
    return render_template("index.html", data_cards=default_data_cards)


@app.route('/download/<filename>')
def download_file(filename):
    return send_file(os.path.join(RESULTS_FOLDER, filename), as_attachment=True)


@app.route('/login', methods=['POST'])
def login():
    username = request.form.get('username')
    password = request.form.get('password')
    user = User.query.filter_by(username=username).first()  # busca el usuario en la base de datos

    if user and user.check_password(password):  # si el usuario existe y la contraseña es correcta
        login_user(user)
        return redirect(url_for('home'))
    else:
        return render_template(
            'index.html', 
            error="Invalid username or password",
            show_login_modal=True
            )

@app.route('/logout', methods=['POST'])
@login_required
def logout():
    logout_user()
    return redirect(url_for('home'))

@app.route('/register', methods=['POST'])
def register():
    username = request.form.get('username')
    password = request.form.get('password')
    confirm_password = request.form.get('confirm_password')
    email = request.form.get('email')

    existing_user = User.query.filter_by(username=username).first()
    if existing_user:
        return render_template(
            'index.html', 
            error="Username already exists",
            show_register_modal=True
            )
        
    if len(password) < 8:
        return render_template(
            'index.html',
            error="Password must be at least 8 characters long.",
            show_register_modal=True
        )

    if password != confirm_password:
        return render_template(
            'index.html',
            error="Passwords do not match.",
            show_register_modal=True
            )
        
    if User.query.filter_by(email=email).first():
        return render_template(
            'index.html', 
            error="Email already registered", 
            show_register_modal=True
            )

    new_user = User(username=username, email=email)
    new_user.set_password(password)  # hashes and stores the password
    db.session.add(new_user)
    db.session.commit()

    login_user(new_user)
    flash("User registered successfully. Welcome!", "success")
    return redirect(url_for('home'))


if __name__ == '__main__':
    with app.app_context():
        db.create_all()
    app.run(debug=True) 