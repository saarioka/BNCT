from io import StringIO
import pandas as pd

def read_wcsv(filename: str) -> pd.DataFrame:
    """
    Parse custom wcsv::ntuple format into a pandas DataFrame.
    file_content: str (whole file as string)
    """
    with open(filename, 'r') as file:
        file_content = file.read()
    
    lines = file_content.strip().splitlines()

    # Defaults
    separator = ","
    vector_separator = ";"
    columns = []

    # Parse header
    data_lines = []
    for line in lines:
        if line.startswith("#separator"):
            sep_ascii = int(line.split()[1])
            separator = chr(sep_ascii)
        elif line.startswith("#vector_separator"):
            vec_ascii = int(line.split()[1])
            vector_separator = chr(vec_ascii)
        elif line.startswith("#column"):
            # Example: "#column double E"
            col_name = line.split()[-1]
            columns.append(col_name)
        elif not line.startswith("#"):
            data_lines.append(line)

    # Join only data lines
    data_str = "\n".join(data_lines)

    # Read into DataFrame
    df = pd.read_csv(
        StringIO(data_str),
        sep=separator,
        names=columns,
        engine="python"
    )
    return df
