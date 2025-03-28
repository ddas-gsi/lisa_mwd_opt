
def load_mwd_settings(file_path):
    params = {}  # Dictionary to store parameters

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()  # Remove spaces and newlines
            if not line or line.startswith("#"):  # Ignore empty lines and comments
                continue

            # Split key and value using '=' and remove inline comments ('>>')
            parts = line.split(">>")[0].split("=")
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip().strip('"')  # Remove double quotes

                # Convert numbers to int or float if applicable
                if value.replace(".", "", 1).isdigit():  # Handle both int and float
                    value = float(value) if "." in value else int(value)

                params[key] = value  # Store in dictionary

    return params

# Usage Example
# file_path = "./set_mwd.txt"  # Path to your parameter file
# parameters = load_mwd_settings(file_path)

# Access parameters
# print(parameters)  # Print all loaded parameters
# print("Input File:", parameters.get("INPUT_FILE", "Not Found"))
# print("Histogram File Path:", parameters.get("HISTOGRAM_FILE_PATH", "Not Found"))
# print("Channel ID:", parameters.get("channelID", "Not Found"))
# print("MWD Length Start:", parameters.get("MWD_Length_Range.start", "Not Found"))