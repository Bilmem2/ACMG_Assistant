# Test alias functionality
from utils.input_handler import InputHandler

# Test choice with aliases
choices = ['heterozygous', 'homozygous', 'hemizygous', 'unknown']
aliases = {'heterozygous': 'het', 'homozygous': 'hom', 'hemizygous': 'hemi'}

print('Testing alias matching...')
print('Choices:', choices)
print('Aliases:', aliases)

# Test alias matching logic
test_inputs = ['het', 'hom', 'hemi', 'heterozygous', 'HET', 'Het']

for test_input in test_inputs:
    value = test_input.lower()
    
    # Check direct matches first
    direct_match = None
    if value in [choice.lower() for choice in choices]:
        direct_match = next(choice for choice in choices if choice.lower() == value)
    
    # Check aliases
    alias_match = None
    for choice, alias in aliases.items():
        if value == alias.lower():
            alias_match = choice
            break
    
    result = direct_match or alias_match
    print(f'Input: "{test_input}" -> Result: {result}')
