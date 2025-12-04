import subprocess
import os

def test_recommend_annotation():
    model_path = '/Users/luna/Desktop/CRBM/AMAS_proj/Models/Test/190_noanno.xml'
    output_csv = 'recommendations.csv'
    
    # Construct the command to run the recommend_annotation script
    command = [
        'python', 'AMAS/recommend_annotation.py', model_path,
        '--cutoff', '0.6', '--save', 'csv', '--outfile', output_csv
    ]
    
    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)
    
    # Check if the command was successful
    assert result.returncode == 0, f"Error: {result.stderr}"
    
    # Check if the output file was created
    assert os.path.exists(output_csv), "Output file was not created."
    
    # Optionally, you can read the output file and check its contents
    with open(output_csv, 'r') as f:
        content = f.readlines()
        assert len(content) > 0, "Output file is empty."

    print("Test passed successfully.")

# Run the test
if __name__ == '__main__':
    test_recommend_annotation() 