import os
import pandas as pd

def extract_alignment_data(file_path):
    # Store the relevant data, initialize with default values
    data = {
        "Total Reads": 0,
        "Aligned Concordantly 0 Times": 0,
        "Aligned Concordantly Exactly 1 Time": "0 (0.00%)",
        "Aligned Concordantly >1 Times": "0 (0.00%)",
        "Aligned Discordantly 1 Time": "0 (0.00%)",
        "Aligned 0 Times Concordantly or Discordantly": 0,
        "Total Mates": 0,
        "Aligned Mates 0 Times": "0 (0.00%)",
        "Aligned Mates Exactly 1 Time": "0 (0.00%)",
        "Aligned Mates >1 Times": "0 (0.00%)",
        "Overall Alignment Rate (%)": "0.0%"
    }
    
    with open(file_path, 'r') as file:
        for line in file:
            if "reads; of these" in line:
                data["Total Reads"] = int(line.split()[0])  # Total Reads
            elif "aligned concordantly 0 times" in line and ">" not in line:
                data["Aligned Concordantly 0 Times"] = int(line.split()[0])  
            elif "aligned concordantly exactly 1 time" in line:
                count = line.split()[0]
                percentage = extract_percentage(line)
                data["Aligned Concordantly Exactly 1 Time"] = f"{count} ({percentage})"
            elif "aligned concordantly >1 times" in line:
                count = line.split()[0]
                percentage = extract_percentage(line)
                data["Aligned Concordantly >1 Times"] = f"{count} ({percentage})"
            elif "aligned discordantly 1 time" in line:
                count = line.split()[0]
                percentage = extract_percentage(line)
                data["Aligned Discordantly 1 Time"] = f"{count} ({percentage})"
            elif "aligned 0 times concordantly or discordantly" in line:
                data["Aligned 0 Times Concordantly or Discordantly"] = int(line.split()[0])
            elif "mates make up the pairs; of these:" in line:
                data["Total Mates"] = int(line.split()[0])
            elif "aligned 0 times" in line:
                count = line.split()[0]
                percentage = extract_percentage(line)
                data["Aligned Mates 0 Times"] = f"{count} ({percentage})"
            elif "aligned exactly 1 time" in line:
                count = line.split()[0]
                percentage = extract_percentage(line)
                data["Aligned Mates Exactly 1 Time"] = f"{count} ({percentage})"
            elif "aligned >1 times" in line:
                count = line.split()[0]
                percentage = extract_percentage(line)
                data["Aligned Mates >1 Times"] = f"{count} ({percentage})"
            elif "overall alignment rate" in line:
                # Extract the overall alignment rate percentage
                percentage = line.split()[0].strip()  # Get the overall alignment rate
                # Ensure only one '%' is appended
                data["Overall Alignment Rate (%)"] = f"{percentage}%" if not percentage.endswith('%') else percentage

    return data

def extract_percentage(line):
    """Helper function to extract the percentage from a given line."""
    try:
        # Check if there are parentheses in the line
        if '(' in line and ')' in line:
            return line.split('(')[1].split('%')[0].strip() + "%"
        else:
            return "0.00%"  # Default if no percentage found
    except IndexError:
        return "0.00%"  # Default if there's an error in extraction

def process_folder(folder_path, output_path):
    all_data = []
    
    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):  # Assuming the files are .txt
            file_path = os.path.join(folder_path, filename)
            file_data = extract_alignment_data(file_path)
            # Remove the .txt extension and add filename
            file_data["Sample"] = os.path.splitext(filename)[0]
            all_data.append(file_data)
    
    # Create a DataFrame
    df = pd.DataFrame(all_data)
    
    # Move the "Sample" column to the first position
    df = df[["Sample"] + [col for col in df.columns if col != "Sample"]]

    # Save to Excel
    df.to_excel(output_path, index=False)

    # Save to a tab-separated .txt file
    #df.to_csv(output_path, sep='\t', index=False)

# usage
folder_path = r'D:/Documents/RNAseq_projects/alignment_summary'  # Replace with your folder path
output_path = r'D:/Documents/RNAseq_projects/alignment_results.xlsx'  # Replace with your desired output path and format

process_folder(folder_path, output_path)

#Completed
