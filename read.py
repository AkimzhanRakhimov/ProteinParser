from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from io import StringIO  # Импортируем StringIO для работы с строками как с файлами

# Чтение аминокислотных цепей из файла (FASTA)
def read_proteins_from_file(file_path):
    proteins = []
    for record in SeqIO.parse(file_path, "fasta"):
        proteins.append(str(record.seq))  # Сохраняем аминокислотные цепи в список
    return proteins

# Удаление стоп-кодонов из последовательности
def remove_stop_codons(protein_sequence):
    return protein_sequence.replace("*", "")  # Убираем символ * (стоп-кодон)

# Отправка запроса BLAST для каждой цепи
def blast_protein_sequence(protein_sequence):
  
    # Отправляем запрос BLAST для текущей цепи аминокислот
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)  # Используем базу данных "nr" (все белки)
    blast_results = result_handle.read()  # Читаем результат поиска
    return blast_results

# Парсинг результатов BLAST и вывод информации
def parse_blast_results(blast_results):
    # Используем StringIO, чтобы обработать строку как файл
    blast_record = NCBIXML.read(StringIO(blast_results))  # Используем StringIO для работы с результатами как с файлом
    for alignment in blast_record.alignments:
        print(f"\nНайден белок: {alignment.title}")
        for hsp in alignment.hsps:
            print(f"  Идентичность: {hsp.identities} из {hsp.align_length}")
            print(f"  Схожесть (e-value): {hsp.expect}")

# Основная функция
def main():
    # Шаг 1: Чтение аминокислотных цепей из файла
    proteins = read_proteins_from_file("proteins_filtered.txt")  # Замените на путь к вашему файлу с белками

    # Шаг 2: Запуск BLAST для каждой цепи и вывод результатов
    for protein in proteins:
        protein = remove_stop_codons(protein)  # Убираем стоп-кодоны
        print(f"\nЗапуск BLAST для белка: {protein[:50]}...")  # Покажем первые 50 аминокислот для удобства
        blast_results = blast_protein_sequence(protein)  # Отправляем запрос BLAST
        parse_blast_results(blast_results)  # Разбираем и выводим результаты

if __name__ == "__main__":
    main()
