import os
import json

def convert_notebook_to_python(ipynb_path):
    # Define the new file name
    py_path = ipynb_path.replace('.ipynb', '.py')
    
    with open(ipynb_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
        
    py_lines = []
    py_lines.append(f"# {'='*60}\n")
    py_lines.append(f"# Auto-generated from notebook: {os.path.basename(ipynb_path)}\n")
    py_lines.append(f"# {'='*60}\n\n")
    
    for cell in nb.get('cells', []):
        cell_type = cell.get('cell_type')
        source = cell.get('source', [])
        
        if not source:
            continue
            
        if cell_type == 'markdown':
            py_lines.append("# " + "-"*40 + "\n")
            py_lines.append("# MARKDOWN CELL\n")
            py_lines.append("# " + "-"*40 + "\n")
            for line in source:
                # Prefix each markdown line with a comment hash
                py_lines.append(f"# {line}")
            py_lines.append("\n\n")
            
        elif cell_type == 'code':
            py_lines.append("# " + "-"*40 + "\n")
            py_lines.append("# CODE CELL\n")
            py_lines.append("# " + "-"*40 + "\n")
            for line in source:
                py_lines.append(line)
            py_lines.append("\n\n")

    with open(py_path, 'w', encoding='utf-8') as f:
        f.writelines(py_lines)
        
    print(f"Successfully converted: {os.path.basename(ipynb_path)} -> {os.path.basename(py_path)}")

if __name__ == "__main__":
    print("Starting notebook conversion...\n")
    
    # Walk through all directories and subdirectories
    for root, dirs, files in os.walk('.'):
        # Skip the archive folder we made earlier so we don't process dead code!
        if 'notebooks_archive' in root or '.ipynb_checkpoints' in root:
            continue
            
        for file in files:
            if file.endswith('.ipynb'):
                full_path = os.path.join(root, file)
                convert_notebook_to_python(full_path)
                
    print("\nAll done!")