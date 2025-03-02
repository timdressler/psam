from psychopy import gui
import os
import random
import pandas as pd
import numpy as np

# Function to ask for Subject ID using PsychoPy's dialog box
def ask_subject_id():
    dlg = gui.Dlg(title="Subject ID Input")
    dlg.addField("Enter Subject ID:")
    dlg.show()
    
    if dlg.OK:
        subject_id = dlg.data[0]
        subject_id = f"sub-{subject_id.zfill(3)}"
        return subject_id
    else:
        return None

# Ask for Subject ID
subject_id = ask_subject_id()

if subject_id is None:
    print("Subject ID input was canceled. Exiting script.")
    exit()

# Set up path for stimuli
script_path = os.path.abspath(__file__)
main_path = os.path.dirname(os.path.dirname(script_path))  # Adjust for actual script location
stimuli_path = os.path.join(main_path, 'data', subject_id, 'stimuli')

# Ensure the directory exists
os.makedirs(stimuli_path, exist_ok=True)

# Maximum consecutive occurrences of the same trial type
max_consecutive = 4  

def create_miniblock():
    """Creates a miniblock with 8 Sham and 8 Probe trials."""
    # Probe trials combinations (Active and Passive)
    active_combinations = [
        ('Active', 'Yes', 'Early', 'Normal'),
        ('Active', 'Yes', 'Early', 'Pitch'),
        ('Active', 'Yes', 'Late', 'Normal'),
        ('Active', 'Yes', 'Late', 'Pitch')
    ]
    
    passive_combinations = [
        ('Passive', 'Yes', 'Early', 'Normal'),
        ('Passive', 'Yes', 'Early', 'Pitch'),
        ('Passive', 'Yes', 'Late', 'Normal'),
        ('Passive', 'Yes', 'Late', 'Pitch')
    ]
    
    # Sham trials (no Probe, no Task)
    sham_trials = ['Sham'] * 8
    
    # Create the miniblock by randomly assigning combinations to active and passive Probe trials
    random.shuffle(active_combinations)
    random.shuffle(passive_combinations)
    
    # The miniblock will contain a shuffled order of Sham and Probe trials
    trials = sham_trials + active_combinations + passive_combinations
    
    # Shuffle to randomize order of the Sham and Probe trials (ensuring no more than 4 consecutive same trials)
    while True:
        random.shuffle(trials)
        if all(trials[i] != trials[i + max_consecutive] for i in range(len(trials) - max_consecutive)):
            break
    
    return trials

def assign_active_passive(trials):
    """Assigns 50% Active and 50% Passive to both Sham and Probe trials, also adding SubjectID column and new probe-related columns."""
    # Start with creating a DataFrame for trials
    df = pd.DataFrame({'Probe': trials})  # Column for Probe or Sham trials
    
    # Add the SubjectID column as the first column
    df.insert(0, 'subj', subject_id)  # Insert the SubjectID as the first column
    
    # Assign the 'Task' column to Active or Passive for the Probe trials
    df['Task'] = np.nan
    df.loc[df['Probe'] == 'Sham', 'Task'] = 'Sham'
    df.loc[df['Probe'] == 'Yes', 'Task'] = random.choices(['Active', 'Passive'], k=len(df[df['Probe'] == 'Yes']))
    
    # Assign probe_type and probe_onset_cat columns to the Probe trials (not Sham)
    df['probe_type'] = np.nan
    df['probe_onset_cat'] = np.nan
    
    # Now, assign the specific combinations (Active/Passive, Early/Late, Normal/Pitch)
    for i, row in df[df['Probe'] == 'Yes'].iterrows():
        if row['Task'] == 'Active':
            # Assign the 4 active combinations randomly
            if row.name < 4:  # First 4 Active probes
                df.at[i, 'probe_type'] = random.choice(['Normal', 'Pitch'])
                df.at[i, 'probe_onset_cat'] = random.choice(['Early', 'Late'])
            else:
                df.at[i, 'probe_type'] = random.choice(['Normal', 'Pitch'])
                df.at[i, 'probe_onset_cat'] = random.choice(['Early', 'Late'])
        elif row['Task'] == 'Passive':
            # Assign the 4 passive combinations randomly
            if row.name < 4:  # First 4 Passive probes
                df.at[i, 'probe_type'] = random.choice(['Normal', 'Pitch'])
                df.at[i, 'probe_onset_cat'] = random.choice(['Early', 'Late'])
            else:
                df.at[i, 'probe_type'] = random.choice(['Normal', 'Pitch'])
                df.at[i, 'probe_onset_cat'] = random.choice(['Early', 'Late'])
    
    return df

def generate_trials(num_trials):
    """Generates and saves the full trial sequence."""
    if num_trials % 16 != 0:
        raise ValueError("Number of trials must be divisible by 16.")
    
    all_trials = []
    
    for _ in range(num_trials // 16):
        miniblock = create_miniblock()
        df_miniblock = assign_active_passive(miniblock)
        all_trials.append(df_miniblock)
    
    final_df = pd.concat(all_trials, ignore_index=True)
    
    conditions_filename = os.path.join(stimuli_path, f"{subject_id}_conditions_shuffled.xlsx")
    final_df.to_excel(conditions_filename, index=False)
    print(f"Trial sequence saved to {conditions_filename}")

# Define the number of trials (modifiable, must be divisible by 16)
num_trials = 64  
generate_trials(num_trials)
print("Test complete")
