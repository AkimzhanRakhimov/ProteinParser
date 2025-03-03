import requests
import json
import re
import os
from datetime import datetime
import tkinter as tk
from tkinter import ttk, messagebox

def extract_numbers(text):
    """Ищет числа в тексте и возвращает список значений"""
    numbers = re.findall(r'\d+\.\d+|\d+', text)
    return [float(n) for n in numbers]

def get_protein_data(uniprot_id):
    url = f"http://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
    except requests.RequestException as e:
        print(f"Ошибка при запросе UniProt ({uniprot_id}): {e}")
        return None

    # Функция белка
    function = next(
        (text["value"] for comment in data.get("comments", [])
         if comment.get("commentType") == "FUNCTION"
         for text in comment.get("texts", [])), 
        "Unknown"
    )

    # Кофакторы
    cofactors = [
        text["value"] for comment in data.get("comments", [])
        if comment.get("commentType") == "PTM"
        for text in comment.get("texts", [])
    ]

    # Локализация
    localization = [
        loc.get("location", {}).get("value", "Unknown") for comment in data.get("comments", [])
        if comment.get("commentType") == "SUBCELLULAR LOCATION"
        for loc in comment.get("subcellularLocations", [])
    ]

    # Взаимодействия
    interactions = [
        interact.get("interactantTwo", {}).get("uniProtKBAccession", "Unknown")
        for comment in data.get("comments", [])
        if comment.get("commentType") == "INTERACTION"
        for interact in comment.get("interactions", [])
    ]

    # Регуляция
    regulation = [
        text["value"] for comment in data.get("comments", [])
        if comment.get("commentType") in ["ACTIVITY REGULATION", "INDUCTION"]
        for text in comment.get("texts", [])
    ]

    # EC Number
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

    # Масса белка
    mass = data.get("sequence", {}).get("mass", "Unknown")
    if mass == "Unknown":
        seq_length = data.get("sequence", {}).get("length", None)
        if seq_length:
            mass = round(seq_length * 110, 2)  # Средняя масса АК ≈ 110 Да

    protein_info = {
        "ID": uniprot_id,
        "Function": function,
        "Cofactors": cofactors,
        "Localization": localization,
        "Interactions": interactions,
        "Regulation": regulation,
        "EC Number": ec_numbers,
        "Mass": mass,
        "Lifetime": get_lifetime(mass),
        "Diffusion": calculate_diffusion(mass),
    }

    # Сохранение в файл
    save_protein_data(protein_info)

    return protein_info

def save_protein_data(protein_info):
    """Сохраняет данные в файл с уникальным именем по дате."""
    folder = "protein_data"
    os.makedirs(folder, exist_ok=True)  # Создаёт папку, если её нет

    filename = datetime.now().strftime("%Y-%m-%d_%H-%M-%S.json")
    filepath = os.path.join(folder, filename)

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(protein_info, f, indent=4, ensure_ascii=False)

    print(f"✅ Данные сохранены в {filepath}")

def get_lifetime(mass):
    """Эмпирическая формула для времени жизни белка."""
    if mass == "Unknown" or not isinstance(mass, (int, float)):
        return "Unknown"
    return round(10 ** (5 - (mass / 50000)), 2)

def calculate_diffusion(mass):
    """Расчёт скорости диффузии по уравнению Стокса-Эйнштейна."""
    if mass == "Unknown" or not isinstance(mass, (int, float)):
        return "Unknown"
    return round(1 / (mass ** 0.5), 5)

# === UI на Tkinter ===
def on_search():
    """Запуск парсинга при вводе UniProt ID."""
    uniprot_id = entry_id.get().strip()
    if not uniprot_id:
        messagebox.showerror("Ошибка", "Введите UniProt ID!")
        return
    
    protein_data = get_protein_data(uniprot_id)
    if not protein_data:
        messagebox.showerror("Ошибка", "Не удалось получить данные!")
        return

    # Очищаем таблицу
    for row in table.get_children():
        table.delete(row)

    # Заполняем таблицу
    for key, value in protein_data.items():
        table.insert("", "end", values=(key, str(value)))

    messagebox.showinfo("Успех", "Данные успешно загружены и сохранены!")

# Создание UI
root = tk.Tk()
root.title("Protein Parser")

# Ввод UniProt ID
frame = tk.Frame(root)
frame.pack(pady=10)
tk.Label(frame, text="UniProt ID:").pack(side="left")
entry_id = tk.Entry(frame, width=20)
entry_id.pack(side="left", padx=5)
tk.Button(frame, text="Парсить", command=on_search).pack(side="left")

# Таблица
columns = ("Параметр", "Значение")
table = ttk.Treeview(root, columns=columns, show="headings")
table.heading("Параметр", text="Параметр")
table.heading("Значение", text="Значение")
table.pack(expand=True, fill="both", padx=10, pady=10)

root.mainloop()
