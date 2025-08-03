from Bio import SeqIO

def analyze_fasta(file_path):
    # Read the first (or only) record from the FASTA file
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

    rev_comp = str(seq.reverse_complement())
    rna = str(seq.transcribe())
    try:
        protein = str(seq.translate(to_stop=True))
    except Exception as e:
        protein = f"Translation error: {e}"

    return {
        "A count": a_count,
        "C count": c_count,
        "T count": t_count,
        "G count": g_count,
        "Length": length,
        "AT/GC Ratio": at_gc_ratio,
        "GC Content %": gc_content,
        "Reverse Complement": rev_comp,
        "Transcription (RNA)": rna,
        "Translation (Protein)": protein
    }


result = analyze_fasta(r"C:\Users\Usuario\Downloads\U00096.3.fasta")
for key, value in result.items():
    print(f"{key}: {value}")
