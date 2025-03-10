Project Title: Protein Data Parser

**Description:**
This project is a Python-based application that retrieves, processes, and stores protein data from the UniProt database. The program accepts a UniProt ID, retrieves relevant protein data (such as function, cofactors, interactions, etc.), and calculates additional properties like diffusion rate and protein lifetime based on empirical formulas. The data is then saved in a JSON file for further analysis. 

**Key Features:**
1. **Data Retrieval**: Fetches protein data from the UniProt database using their REST API in JSON format.
2. **Data Extraction**: Extracts key protein-related information such as function, cofactors, subcellular localization, interactions, regulation, EC numbers, and mass.
3. **Advanced Calculations**: Computes protein lifetime using an empirical formula and diffusion rate using the Stokes-Einstein equation, based on the protein mass.
4. **Data Saving**: Automatically saves the extracted protein data in a uniquely timestamped JSON file.
5. **User Interface**: Built using Tkinter for a simple graphical interface that allows users to input a UniProt ID, display the retrieved data, and save it for future use.

**Technologies Used:**
- Python (requests, re, json, os, datetime, Tkinter)
- API integration (UniProt REST API)
- Data processing and visualization (JSON manipulation)

**Skills Demonstrated:**
- **API Integration**: Used the UniProt API to retrieve detailed protein information.
- **Data Parsing and Processing**: Extracted and processed data from JSON responses and performed calculations based on protein properties.
- **User Interface Design**: Designed a simple graphical user interface (GUI) for ease of use using Tkinter.
- **Error Handling**: Implemented error handling for potential issues during data retrieval and processing.

**How to Run:**
1. Install Python 3.x.
2. Install required libraries using:
   ```
   pip install requests
   ```
3. Run the script `protein_parser.py`.
4. Enter a UniProt ID when prompted to fetch the protein data(For example:P000742).

