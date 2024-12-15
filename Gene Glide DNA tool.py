import tkinter as tk
from tkinter import filedialog, messagebox  # Importing filedialog and messagebox

# Codon Table for DNA to Protein Translation
codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
    'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCT': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': '', 'UAG': '',
    'UGC': 'C', 'UGU': 'C', 'UGA': '', 'UGG': 'W',
}

def translate_sequence(sequence, is_rna=False):
    """Translate a DNA or RNA sequence into its corresponding amino acid sequence."""
    amino_acids = ""
    if is_rna:
        sequence = sequence.replace('U', 'T')  # Convert RNA to DNA for translation

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        amino_acid = codon_table.get(codon, '_')
        amino_acids += amino_acid

    return amino_acids

def translate_file(file_path):
    """Translate the sequence from a file."""
    try:
        with open(file_path, 'r') as file:
            sequence = file.read().strip()
            return translate_sequence(sequence.upper())
    except FileNotFoundError:
        return "Error: File not found."
    except Exception as e:
        return f"Error: {str(e)}"

def browse_file():
    """Open a file dialog to select a FASTA file."""
    try:
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", ".fasta;.fa")])
        if file_path:
            with open(file_path, 'r') as file:
                # Read the FASTA file and extract the sequence
                lines = file.readlines()
                sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))  # Exclude FASTA headers
                entry_text.delete("1.0", tk.END)  # Clear the current content
                entry_text.insert(tk.END, sequence)  # Insert the sequence into the input text box
    except Exception as e:
        error_label.config(text=f"Error: {str(e)}", fg="red")

def button1_click():
    """Handle the click event for translating the sequence."""
    input_seq = entry_text.get("1.0", "end-1c").upper()
    is_rna = dna_var.get() == 1  # 1 indicates RNA

    if not is_rna and 'U' in input_seq:
        error_label.config(text="Error: Enter a valid DNA sequence.", fg="red")
        return
    if is_rna and 'T' in input_seq:
        error_label.config(text="Error: Enter a valid RNA sequence.", fg="red")
        return
        
    try:
        if is_rna:
            # Translate RNA sequence to protein
            amino_acid_sequence = translate_sequence(input_seq, is_rna=True)
            sequence_text.config(state=tk.NORMAL)
            sequence_text.delete("1.0", tk.END)
            sequence_text.insert(tk.END, "Amino acid sequence (RNA to Protein): ", "bold_red")
            sequence_text.insert(tk.END, f"{amino_acid_sequence}\n")
        else:
            amino_acid_sequence = translate_sequence(input_seq, is_rna=False)
            sequence_text.config(state=tk.NORMAL)
            sequence_text.delete("1.0", tk.END)
            sequence_text.insert(tk.END, "Amino acid sequence (DNA to Protein): ", "bold_blue")
            sequence_text.insert(tk.END, f"{amino_acid_sequence}\n")

        sequence_text.config(state=tk.DISABLED)
        error_label.config(text="", fg="red")  # Clear error message if successful
    except Exception as e:
        error_label.config(text=f"Error: {str(e)}", fg="red")

def save_result():
    """Save the translation result to a file."""
    try:
        result = sequence_text.get("1.0", tk.END).strip()  # Get the result text
        if not result:
            messagebox.showwarning("Warning", "There is no result to save.")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        if file_path:
            with open(file_path, 'w') as file:
                file.write(result)
            messagebox.showinfo("Success", "Result saved successfully.")
    except Exception as e:
        error_label.config(text=f"Error: {str(e)}", fg="red")

def reset_fields():
    """Reset input and output fields."""
    entry_text.delete("1.0", tk.END)  # Clear input text
    sequence_text.config(state=tk.NORMAL)
    sequence_text.delete("1.0", tk.END)  # Clear output text
    sequence_text.config(state=tk.DISABLED)  # Disable output text box
    error_label.config(text="", fg="red")  # Clear error message
    dna_var.set(0)  # Reset the radio button selection

# Create the main application window
root = tk.Tk()
root.title("GeneGlide")
root.geometry("800x600")
root.configure(bg='#f0f0f0')  # Light background color

# Header Section
header_frame = tk.Frame(root, bg='#3b5998')
header_frame.pack(fill=tk.X)

# Tool Name and Icon
title_label = tk.Label(header_frame, text="GeneGlide", font=("Morro", 36, "bold"), bg='#3b5998', fg="white")
title_label.pack(side=tk.LEFT, padx=10, pady=10)

# Add DNA Icon (using text here, you can replace with an actual icon)
dna_icon = tk.Label(header_frame, text="🧬", font=("Arial", 36), bg='#3b5998', fg="white")
dna_icon.pack(side=tk.LEFT)

# Ribbon Section for Buttons
ribbon_frame = tk.Frame(header_frame, bg='#3b5998')
ribbon_frame.pack(side=tk.RIGHT)

def open_features():
    """Open App Features in a new window."""
    features_window = tk.Toplevel(root)
    features_window.title("App Features")
    features_content = "This tool translates DNA/RNA sequences to amino acids.\n\nFeatures:\n- Input DNA or RNA sequences\n- Browse FASTA files\n- Save translation results"
    label = tk.Label(features_window, text=features_content, padx=10, pady=10)
    label.pack()

def open_guide():
    """Open Tool Guide in a new window."""
    guide_window = tk.Toplevel(root)
    guide_window.title("Tool Guide")
    guide_content = "To use this tool:\n1. Enter a DNA or RNA sequence in the input box.\n2. Select DNA or RNA.\n3. Click 'Translate' to see the amino acid sequence."
    label = tk.Label(guide_window, text=guide_content, padx=10, pady=10)
    label.pack()

def contact_developer():
    """Open Contact Us information in a new window."""
    contact_window = tk.Toplevel(root)
    contact_window.title("Contact Developer")
    contact_content = "Contact Developer:\n\nNoor SKhan.\nEmail: emankhan7228@gmail.com"
    label = tk.Label(contact_window, text=contact_content, padx=10, pady=10)
    label.pack()

features_button = tk.Button(ribbon_frame, text="App Features", command=open_features)
features_button.pack(side=tk.LEFT, padx=5, pady=5)

guide_button = tk.Button(ribbon_frame, text="Tool Guide", command=open_guide)
guide_button.pack(side=tk.LEFT, padx=5, pady=5)

contact_button = tk.Button(ribbon_frame, text="Contact Developer", command=contact_developer)
contact_button.pack(side=tk.LEFT, padx=5, pady=5)

# Input Section
input_frame = tk.Frame(root, bg='#f0f0f0')
input_frame.pack(padx=10, pady=10)

input_label = tk.Label(input_frame, text="Enter DNA or RNA Sequence:", bg='#f0f0f0')
input_label.pack()

entry_text = tk.Text(input_frame, height=4, width=60)
entry_text.pack(pady=5)

# Option to select input type (DNA/RNA)
dna_var = tk.IntVar()
dna_radio = tk.Radiobutton(input_frame, text="DNA", variable=dna_var, value=0, bg='#f0f0f0')
dna_radio.pack(side=tk.LEFT, padx=5)

rna_radio = tk.Radiobutton(input_frame, text="RNA", variable=dna_var, value=1, bg='#f0f0f0')
rna_radio.pack(side=tk.LEFT, padx=5)

browse_button = tk.Button(input_frame, text="Browse FASTA File", command=browse_file)
browse_button.pack(side=tk.LEFT, padx=5)

# Button Section
button_frame = tk.Frame(root, bg='#f0f0f0')
button_frame.pack(pady=10)

translate_button = tk.Button(button_frame, text="Translate", command=button1_click)
translate_button.pack(side=tk.LEFT, padx=5)

save_button = tk.Button(button_frame, text="Save Result", command=save_result)
save_button.pack(side=tk.LEFT, padx=5)

reset_button = tk.Button(button_frame, text="Reset", command=reset_fields)
reset_button.pack(side=tk.LEFT, padx=5)

# Output Section
output_frame = tk.Frame(root, bg='#f0f0f0')
output_frame.pack(padx=10, pady=10)

output_label = tk.Label(output_frame, text="Amino Acid Sequence:", bg='#f0f0f0')
output_label.pack()

sequence_text = tk.Text(output_frame, height=10, width=60, state=tk.DISABLED)
sequence_text.pack(pady=5)

# Error message display
error_label = tk.Label(root, text="", bg='#f0f0f0', fg="red")
error_label.pack(pady=5)

# Start the main event loop
root.mainloop()
