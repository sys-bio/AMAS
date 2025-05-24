import json
import re

test = """
RESULT:
{"CaI": ["calcium ion", "Ca2+", "calcium(2+)"], 
"IP3": ["inositol 1,4,5-trisphosphate", "D-myo-inositol 1,4,5-trisphosphate", "IP3(4-)"], 
"CaO": ["UNK"], 
"CaS": ["calcium ion", "Ca2+", "calcium(2+)"], 
"CaM": ["calmodulin", "calcium-modulated protein", "CaM(1-)]}
Reason: Based on the model name "Abell2011_CalciumSignaling_WithoutAdaptation" and the reactions, it is likely that CaI and CaS represent calcium ions in different locations (e.g., cytosol and endoplasmic reticulum). IP3 is a well-known second messenger involved in calcium signaling. CaM is likely calmodulin, a calcium-binding protein. CaO is unknown as there is no clear indication of its identity in the provided information."""

# Test direct JSON parsing
result_match = re.search(r'RESULT:(.*?)Reason:', test, re.DOTALL)
if result_match:
    result_str = result_match.group(1).strip()
    # Find the JSON object in the result
    json_match = re.search(r'({.*})', result_str, re.DOTALL)
    if json_match:
        json_str = json_match.group(1)
        print(f'JSON string: {json_str}')
        try:
            synonyms_dict = json.loads(json_str)
            print(f'Parsed successfully: {synonyms_dict}')
        except json.JSONDecodeError as e:
            print(f'JSON decode error: {e}')
            
            # Count brackets to check if properly balanced
            open_braces = json_str.count('{')
            close_braces = json_str.count('}')
            print(f"Open braces: {open_braces}, Close braces: {close_braces}")
            
            # Fix missing closing bracket if needed
            if open_braces > close_braces:
                fixed_json = json_str + "}"
                print(f"Fixed JSON: {fixed_json}")
                try:
                    synonyms_dict = json.loads(fixed_json)
                    print(f'Fixed version parsed successfully: {synonyms_dict}')
                except json.JSONDecodeError as e:
                    print(f'Fixed JSON still has error: {e}')
            
            # Look for the problematic part - missing quote for CaM entry
            print("\nChecking for specific issues:")
            last_entry = re.search(r'"CaM": \[(.*?)\]\}', json_str)
            if last_entry:
                last_entry_content = last_entry.group(1)
                print(f"Last entry content: {last_entry_content}")
                
                # Check if the last item in the array is missing a closing quote
                if "CaM(1-)" in last_entry_content:
                    print("Found issue: Missing quote after CaM(1-)")
                    
                    # Fix the JSON by adding the missing quote
                    fixed_json = json_str.replace('CaM(1-)]', 'CaM(1-)"]')
                    print(f"Fixed JSON: {fixed_json}")
                    
                    try:
                        synonyms_dict = json.loads(fixed_json)
                        print(f'Fixed version parsed successfully: {len(synonyms_dict)} entries')
                    except json.JSONDecodeError as e:
                        print(f'Fixed JSON still has error: {e}')
            else:
                print("Could not identify the last entry pattern") 