"""
##################################################################################################
# Title: Docking Score Z-Score Analysis Tool
# Developed by: Mine Isaoglu, Ph.D.
# Principal Investigator: Serdar Durdagi, Ph.D.
# Affiliation: Computational Drug Design Center (HITMER), Faculty of Pharmacy,
#              Bahçeşehir University, Istanbul, Turkey
# Version: January 2025
##################################################################################################
"""

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from fpdf import FPDF

def parse_float(value):
    try:
        # If the value contains a comma and has at least two characters after the comma
        if ',' in value and len(value.split(',')[1]) >= 2:
            return float(value.replace(',', '.'))
        # If the value contains a comma but has less than two characters after the comma
        elif ',' in value:
            return float(value.replace(',', '.') + '00')  # Add '00' after the comma
        else:
            return float(value)
    except ValueError:
        return None  # Return None if the value cannot be converted to float

def z_score(data, x):
    mean = np.mean(data)
    std_dev = np.std(data)
    z = (x - mean) / std_dev
    return z

def print_last_molecule_docking(csv_file_path):
    # Open and read the CSV file
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        data = [(row['Title'], row['docking score']) for row in reader]

    # Get the last molecule's docking score
    last_molecule = data[-1]  # Get the last item from the data
    molecule_name, docking_score = last_molecule

    print(f"The docking score of the last molecule ({molecule_name}) is: {docking_score}")
    return float(docking_score)  # Return the score as a float to use it in the plot

def main():
    # Define file names and paths
    file_name = 'file_name.csv'
    updated_file_name = 'file_name_with_Z_Scores.csv'
    output_pdf_name = "Z-Score_Table.pdf"
    output_image_name = "Z_Score_Distribution_Curve.png"
    docking_output_image_name = "Docking_Score_Distribution_Curve.png"

    # Get the working directory
    working_directory = os.getcwd()

    # Update file paths
    file_path = os.path.join(working_directory, file_name)
    updated_file_path = os.path.join(working_directory, updated_file_name)
    output_pdf_path = os.path.join(working_directory, output_pdf_name)
    output_image_path = os.path.join(working_directory, output_image_name)
    docking_output_image_path = os.path.join(working_directory, docking_output_image_name)

    molecule_name_column = 'Title'  # Name of the column containing molecule names
    docking_score_column = 'docking score'  # Name of the column containing docking score values
    z_score_threshold_lower = -1.960  # Z score threshold value (lower bound)

    # Open and read the CSV file
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        data = [(row[molecule_name_column], parse_float(row[docking_score_column])) for row in reader]

    # Filter out None values
    data = [(molecule, value) for molecule, value in data if value is not None]

    # Calculate Z scores
    docking_scores = np.array([item[1] for item in data])
    z_scores = [(x - np.mean(docking_scores)) / np.std(docking_scores) for x in docking_scores]

    # Filter molecules with Z score less than or equal to -1
    selected_data = [(item[0], item[1], '{:.6f}'.format(z)) for item, z in zip(data, z_scores) if z <= z_score_threshold_lower]

    # Update data with Z scores
    data = [(item[0], item[1], z) for item, z in zip(data, z_scores)]

    # Write the updated data to a new CSV file
    with open(updated_file_path, 'w', newline='') as csvfile:
        fieldnames = [molecule_name_column, docking_score_column, 'Z Score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for item in selected_data:  # Only write selected data with Z-score <= -1
            writer.writerow({molecule_name_column: item[0], docking_score_column: item[1], 'Z Score': item[2]})

    # Create PDF
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", style='B', size=12)

    # Write selected molecules' information to PDF
    pdf.cell(200, 10, txt="Compounds with Z-Score <= -1.960", ln=True, align="C")
    pdf.ln(10)  # New line

    # Add total number of selected compounds information
    pdf.cell(200, 10, txt=f"Total Number of Selected Compounds: {len(selected_data)}", ln=True, align="C")
    pdf.ln(10)  # New line

    # Add table headers
    pdf.set_font("Arial", style='B', size=8)
    pdf.cell(60, 10, "Compound ID", border=1)
    pdf.cell(60, 10, "Docking Score (kcal/mol)", border=1)
    pdf.cell(60, 10, "Calculated Z-Score", border=1)
    pdf.ln(10)  # New line

    # Add table data
    pdf.set_font("Arial", size=10)
    for item in selected_data:
        pdf.cell(60, 10, item[0], border=1)
        pdf.cell(60, 10, str(item[1]), border=1)
        pdf.cell(60, 10, item[2], border=1)  # Use the formatted Z score
        pdf.ln(10)  # New line

    # Saving the table output to a file
    pdf.output(output_pdf_path)

    # Plotting the Z-score graph
    plt.figure(figsize=(8, 6))
    plt.title('Z-Score Normal Distribution', fontweight='bold')
    plt.xlabel('Z-Score')
    plt.ylabel('Probability Density')

    # Plotting the normal distribution curve for Z-scores with thicker blue line
    mu, sigma = 0, 1  # Mean and standard deviation
    x = np.linspace(-5, 5, 100)
    plt.plot(x, norm.pdf(x, mu, sigma), color='blue', label='Z-Score Normal Distribution', linewidth=2.5)
    plt.axvline(z_score_threshold_lower, color='red', linestyle='--')  # Adding a vertical line at the lower threshold
    plt.axvline(0, color='red', linestyle='-')  # Adding a vertical line at x=0
    plt.legend()

    # Saving the Z-score graph
    plt.savefig(output_image_path, bbox_inches='tight', dpi=300)
    plt.close()

    # Plotting the docking score distribution graph
    plt.figure(figsize=(8, 6))
    docking_values = [item[1] for item in data]
    plt.title('Docking Score Normal Distribution', fontweight='bold')
    plt.xlabel('Docking Score (kcal/mol)')
    plt.ylabel('Probability Density')

    # Plotting the normal distribution curve for docking scores with thicker green line
    docking_mean = np.mean(docking_values)
    docking_std_dev = np.std(docking_values)
    docking_x = np.linspace(min(docking_values), max(docking_values), 100)
    plt.plot(docking_x, norm.pdf(docking_x, docking_mean, docking_std_dev), color='green', label='Docking Score Normal Distribution', linewidth=2.5)
    plt.legend()

    # Get the exact docking score for the last molecule (from the table)
    last_docking_value = print_last_molecule_docking(updated_file_path)

    # Plotting the exact docking score from the table
    plt.axvline(last_docking_value, color='red', linestyle='--')  # Adding a vertical line at the last molecule's docking score
    plt.text(last_docking_value - 1, 0.01, f'Fitness= {last_docking_value:.2f}', color='red', fontsize=12, fontweight='bold', ha='right')  # Adding text to denote Fitness value with a slight offset

    # Saving the docking score distribution graph
    plt.savefig(docking_output_image_path, bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
