import requests
import json
import re
import os
from datetime import datetime
import tkinter as tk
from tkinter import ttk, messagebox

hydrophobicity_scale = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -0.5, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

def extract_numbers(text):
    """Finds numbers in text and returns a list of values"""
    numbers = re.findall(r'\d+\.\d+|\d+', text)
    return [float(n) for n in numbers]

def calculate_charge(sequence, ph=7.0):
    """Calculates the overall charge of the protein at a given pH"""
    pka_values = {
        'D': 3.9, 'E': 4.1, 'H': 6.0, 'C': 8.3, 'Y': 10.1,
        'K': 10.5, 'R': 12.5, 'N-term': 9.0, 'C-term': 2.0
    }
    
    charge = 0
    if sequence:
        charge += 1 / (1 + 10 ** (ph - pka_values['N-term']))  # N-terminus
        charge -= 1 / (1 + 10 ** (pka_values['C-term'] - ph))  # C-terminus

        for aa in sequence:
            if aa in ['D', 'E', 'C', 'Y']:
                charge -= 1 / (1 + 10 ** (pka_values[aa] - ph))
            elif aa in ['H', 'K', 'R']:
                charge += 1 / (1 + 10 ** (ph - pka_values[aa]))
    
    return round(charge, 2)

def calculate_hydrophobicity(sequence):
    """Calculates the hydrophobicity of the protein based on the sequence"""
    hydrophobicity = 0
    valid_aa_count = 0  # Counting valid amino acids

    if sequence:
        for aa in sequence:
            value = hydrophobicity_scale.get(aa)
            if value is not None:
                hydrophobicity += value
                valid_aa_count += 1

    if valid_aa_count > 0:
        hydrophobicity = hydrophobicity / valid_aa_count  # Average hydrophobicity

    return round(hydrophobicity, 2)

def calculate_radius(mass):
    """Calculates the radius of the protein using the formula: R ≈ 0.066 × (M ^ 0.37)"""
    if mass == "Unknown" or not isinstance(mass, (int, float)):
        return "Unknown"
    radius_nm = 0.066 * (mass ** 0.37)
    return round(radius_nm, 3)  # In nanometers

def predict_shape(length):
    """Predicts the shape of the protein based on its length"""
    if length < 100:
        return "Compact/Spherical"
    elif length < 300:
        return "Medium, possibly with domains"
    else:
        return "Probably elongated/multidomain"

def get_protein_data(uniprot_id, ph):
    url = f"http://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
    except requests.RequestException as e:
        print(f"Error during UniProt request ({uniprot_id}): {e}")
        return None

    # Protein function
    function = next(
        (text["value"] for comment in data.get("comments", [])
         if comment.get("commentType") == "FUNCTION"
         for text in comment.get("texts", [])), 
        "Unknown"
    )
    # Amino acid sequence
    sequence = data.get("sequence", {}).get("value", "Unknown")

    cofactors = [
        text["value"] for comment in data.get("comments", [])
        if comment.get("commentType") == "PTM"
        for text in comment.get("texts", [])
    ]

    # Localization
    localization = [
        loc.get("location", {}).get("value", "Unknown") for comment in data.get("comments", [])
        if comment.get("commentType") == "SUBCELLULAR LOCATION"
        for loc in comment.get("subcellularLocations", [])
    ]

    # Interactions
    interactions = [
        interact.get("interactantTwo", {}).get("uniProtKBAccession", "Unknown")
        for comment in data.get("comments", [])
        if comment.get("commentType") == "INTERACTION"
        for interact in comment.get("interactions", [])
    ]

    # Regulation (activation/inhibition)
    regulation = [
        text["value"] for comment in data.get("comments", [])
        if comment.get("commentType") in ["ACTIVITY REGULATION", "INDUCTION"]
        for text in comment.get("texts", [])
    ]

    # EC Number (searching in multiple places)
    ec_numbers = []
    desc = data.get("proteinDescription", {})

    for section in ["recommendedName", "alternativeNames"]:
        names = desc.get(section, [])
        if not isinstance(names, list):
            names = [names]
        for name in names:
            ec_numbers.extend(ec.get("value") for ec in name.get("ecNumbers", []))

    for comment in data.get("comments", []):
        if comment.get("commentType") == "CATALYTIC ACTIVITY":
            ec = comment.get("reaction", {}).get("ecNumber")
            if ec and ec not in ec_numbers:
                ec_numbers.append(ec)

    # Protein mass (if sequence.mass is missing, calculate by sequence length)
    mass = data.get("sequence", {}).get("mass", "Unknown")
    if mass == "Unknown":
        seq_length = data.get("sequence", {}).get("length", None)
        if seq_length:
            mass = round(seq_length * 110, 2)  # Average AA mass ≈ 110 Da

    # Protein charge
    charge = calculate_charge(sequence, ph) if sequence != "Unknown" else "Unknown"

    hydrophobicity = calculate_hydrophobicity(sequence) if sequence != "Unknown" else "Unknown"
    radius = calculate_radius(mass)
    length = len(sequence) if sequence != "Unknown" else 0
    shape = predict_shape(length)
    diffusion = calculate_diffusion(mass)
    function_result = categorize_function(charge, hydrophobicity, shape, diffusion, cofactors, localization, interactions, regulation, ec_numbers)

    protein_info = {
        "ID": uniprot_id,
        "Sequence": sequence,
        "Mass": mass,
        "Charge at pH 7.0": charge,
        "Hydrophobicity": hydrophobicity,
        "Lifetime": get_lifetime(mass),
        "Radius (nm)": radius,
        "Predicted Shape": shape,
        "Lifetime": get_lifetime(mass),
        "Diffusion": calculate_diffusion(mass),
        "Cofactors": cofactors,
        "Localization": localization,
        "Interactions": interactions,
        "Regulation": regulation,
        "EC Number": ec_numbers,
        "Function": function,
        "Predicted Function": function_result["Predicted Function"],
        "Confidence (%)": function_result["Confidence (%)"],
        "Alternative Functions": function_result["Alternatives"],
        "Prediction Notes": "; ".join(function_result["Comments"]),
    }

    # Save to file
    save_protein_data(protein_info)

    return protein_info

def categorize_function(charge, hydrophobicity, shape, diffusion, cofactor_info, localization, interactions, regulation, EC_number):
    possible_functions = []
    explanation = []

    # Confidence level (0-100)
    confidence = {
        "Enzyme": 0,
        "Signaling": 0,
        "Membrane": 0,
        "Structural": 0,
        "Unknown": 0,
    }

    # === Enzyme ===
    if charge > 4 and hydrophobicity < 0:
        confidence["Enzyme"] += 60
        explanation.append("High positive charge and low hydrophobicity are characteristic of enzymes.")
    if shape == "Compact/Spherical":
        confidence["Enzyme"] += 20
        explanation.append("Compact shape increases the likelihood of enzymatic activity.")
    if EC_number:
        confidence["Enzyme"] += 30
        explanation.append(f"EC number present ({EC_number}) — supports enzymatic activity.")

    # === Signaling ===
    if diffusion > 0.005:
        confidence["Signaling"] += 50
        explanation.append("High diffusion could indicate a signaling role.")
    if shape == "Compact/Spherical" and charge < 2:
        confidence["Signaling"] += 30
        explanation.append("Small neutral protein with a spherical shape — possible signaling peptide.")
    if "Secreted" in localization:
        confidence["Signaling"] += 20
        explanation.append("Protein is secreted, which may indicate a signaling role.")

    # === Membrane ===
    if hydrophobicity > 2.0:
        confidence["Membrane"] += 70
        explanation.append("High hydrophobicity suggests possible membrane integration.")
    if shape != "Compact/Spherical":
        confidence["Membrane"] += 10
    if any("Membrane" in interaction for interaction in interactions):
        confidence["Membrane"] += 30
        explanation.append("Protein interacts with membrane proteins.")

    # === Structural ===
    if shape == "Probably elongated/multidomain":
        confidence["Structural"] += 60
        explanation.append("Long, elongated proteins often serve a structural role.")
    if diffusion < 0.002:
        confidence["Structural"] += 20
        explanation.append("Low diffusion may indicate a structural role.")
    if "Inhibited by" in regulation:
        confidence["Structural"] += 10
        explanation.append("Inhibition regulates protein function, suggesting a structural role.")

    # === Function selection ===
    best_guess = max(confidence, key=confidence.get)
    score = confidence[best_guess]

    # If confidence is too low
    if score < 30:
        best_guess = "Unknown"
        explanation.append("Not enough evidence for reliable function prediction.")

    return {
        "Predicted Function": best_guess,
        "Confidence (%)": score,
        "Alternatives": {k: v for k, v in confidence.items() if k != best_guess and v > 10},
        "Comments": explanation
    }

def save_protein_data(protein_info):
    """Saves protein data to a file with a unique name based on the date."""
    folder = "protein_data"
    os.makedirs(folder, exist_ok=True)  # Creates folder if it doesn't exist

    filename = datetime.now().strftime("%Y-%m-%d_%H-%M-%S.json")
    filepath = os.path.join(folder, filename)

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(protein_info, f, indent=4, ensure_ascii=False)

    print(f"✅ Data saved to {filepath}")

def get_lifetime(mass):
    """Empirical formula for protein lifetime."""
    if mass == "Unknown" or not isinstance(mass, (int, float)):
        return "Unknown"
    return round(10 ** (5 - (mass / 50000)), 2)

def calculate_diffusion(mass):
    """Calculates diffusion rate using the Stokes-Einstein equation."""
    if mass == "Unknown" or not isinstance(mass, (int, float)):
        return "Unknown"
    return round(1 / (mass ** 0.5), 5)

# === UI with Tkinter ===
def on_search():
    """Starts parsing when UniProt ID and pH are entered."""
    uniprot_id = entry_id.get().strip()
    try:
        ph = float(entry_ph.get().strip())  # Gets pH value from the input field
    except ValueError:
        messagebox.showerror("Error", "Enter a valid pH value!")
        return
    
    if not uniprot_id:
        messagebox.showerror("Error", "Enter a UniProt ID!")
        return
    
    protein_data = get_protein_data(uniprot_id, ph)
    if not protein_data:
        messagebox.showerror("Error", "Failed to retrieve data!")
        return

    # Clear the table
    for row in table.get_children():
        table.delete(row)

    # Fill the table
    for key, value in protein_data.items():
        table.insert("", "end", values=(key, value))

root = tk.Tk()
root.title("Protein Information Search")

# UI Elements
frame = ttk.Frame(root, padding="10")
frame.pack(fill="both", expand=True)

# UniProt ID input
label_id = ttk.Label(frame, text="UniProt ID:")
label_id.grid(row=0, column=0, sticky="w", padx=5, pady=5)
entry_id = ttk.Entry(frame, width=20)
entry_id.grid(row=0, column=1, padx=5, pady=5)

# pH input
label_ph = ttk.Label(frame, text="pH (default 7.0):")
label_ph.grid(row=1, column=0, sticky="w", padx=5, pady=5)
entry_ph = ttk.Entry(frame, width=20)
entry_ph.grid(row=1, column=1, padx=5, pady=5)

# Search button
search_button = ttk.Button(frame, text="Search", command=on_search)
search_button.grid(row=2, column=0, columnspan=2, pady=10)

# Results Table
columns = ["Property", "Value"]
table = ttk.Treeview(frame, columns=columns, show="headings", height=10)
table.grid(row=3, column=0, columnspan=2, padx=5, pady=5)

for col in columns:
    table.heading(col, text=col)

root.mainloop()
