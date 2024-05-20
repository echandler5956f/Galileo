import torch
from transformers import GPT2Tokenizer, GPT2Model
from sklearn.preprocessing import OneHotEncoder
from keras.preprocessing.sequence import pad_sequences

# We want to convert optimal control problem inputs (task descriptions, knot times, contact sequences, initial and final states) into a format that can be fed into a neural network.
# These optimal control problem inputs are in the context of legged robot locomotion, where 
#   - task descriptions are high-level commands (e.g., "walk", "run", "jump"), 
#   - knot times are the time intervals for the trajectory optimization, 
#   - contact sequences are binary values indicating contact or no contact with the ground, 
#   - initial and final states are the robot configurations at the start and end of the trajectory

# The goal is to predict the best initial guess for the optimal control problem based on these inputs.

# Example batch of sentences
task_descriptions_batch = ["walk", "jump", "run"]

knot_times_batch = [
    [0.5, 0.4, 0.4, 0.4, 0.4, 0.5],
    [0.2, 0.2, 0.2],
    [0.3, 0.3, 0.3, 0.3, 0.3],
]

# Each element in the contact_sequence_batch is a contact sequences for a single task. Each contact sequence is a list of binary values where 0 indicates contact and -1 indicates no contact.
contact_sequences_batch = [
    [[0, 0, -1, -1], [-1, 0, -1, -1], [0, -1, -1, -1], [-1, 0, -1, -1], [0, -1, -1, -1], [0, 0, -1, -1]],
    [[0, 0, -1, -1], [-1, -1, -1, -1], [0, 0, -1, -1]],
    [[0, 0, -1, -1], [-1, 0, -1, -1], [0, -1, -1, -1], [-1, 0, -1, -1], [0, 0, -1, -1]],
]

q0_batch = [
    [0, 0, 0.864977, 0, 0, 0, 1, 0, 0, 0, -0.5, -1.2, 2, 0.75, 0, 0.25, 0, 0, 0.5, 1.2, 2, -0.75, 0, -0.25, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0],
    [0, 0, 0.864977, 0, 0, 0, 1, 0, 0, 0, 0, -1.5, 0, 0, 0, 0, 0, 0, 0, 1.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0],
    [0, 0, 0.864977, 0, 0, 0, 1, 0, 0, 0, -0.5, -1.2, 2, 0.75, 0, 0.25, 0, 0, 0.5, 1.2, 2, -0.75, 0, -0.25, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0]
]

qf_batch = [
    [0.25, 0, 0.864977, 0, 0, 0, 1, 0, 0, 0, -0.5, -1.2, 2, 0.75, 0, 0.25, 0, 0, 0.5, 1.2, 2, -0.75, 0, -0.25, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0],
    [0.1, 0, 0.864977, 0, 0, 0, 1, 0, 0, 0, 0, -1.5, 0, 0, 0, 0, 0, 0, 0, 1.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0],
    [0.375, 0, 0.864977, 0, 0, 0, 1, 0, 0, 0, -0.5, -1.2, 2, 0.75, 0, 0.25, 0, 0, 0.5, 1.2, 2, -0.75, 0, -0.25, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0]
]

# Initialize the GPT2 tokenizer and model
tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
model = GPT2Model.from_pretrained("gpt2")

# Tokenize the sentences and convert to input IDs
input_ids = [
    tokenizer.encode(text, add_special_tokens=True) for text in task_descriptions_batch
]

# Pad the input IDs so that all sequences are the same length
input_ids = pad_sequences(
    input_ids, maxlen=100, dtype="long", value=0, truncating="post", padding="post"
)

# Convert input IDs to a PyTorch tensor
input_ids = torch.tensor(input_ids)

# Get the GPT2 embeddings
with torch.no_grad():
    embeddings = model(input_ids)

# `embeddings` is a tuple where the first element is the hidden states from the last layer of the model.
# We take the first token (the [CLS] token) from the last layer for each sentence.
sentence_embeddings = [
    embeddings[0][i, 0, :].numpy() for i in range(len(task_descriptions_batch))
]

print(sentence_embeddings)

# Initialize a one-hot encoder
encoder = OneHotEncoder(sparse=False)

# Flatten the contact sequences
flattened_contact_sequences = [sequence for batch in contact_sequences_batch for sequence in batch]

# Fit the encoder to the data and transform the data
one_hot_contact_sequences = encoder.fit_transform(flattened_contact_sequences)

# Determine the number of rows and columns for the reshaping
num_rows = len(contact_sequences_batch)
num_cols = max(len(seq) for batch in contact_sequences_batch for seq in batch)

# Reshape the data
one_hot_contact_sequences = one_hot_contact_sequences.reshape(num_rows, num_cols, -1)

import numpy as np
from torch.nn.utils.rnn import pad_sequence

# Assume the following variables are your processed data
task_descriptions = torch.tensor(sentence_embeddings)  # shape: (batch_size, embedding_size)
knot_times = pad_sequence([torch.tensor(k) for k in knot_times_batch], batch_first=True, padding_value=-1) # shape: (batch_size, max_sequence_length)

# Determine the maximum sequence length
max_sequence_length = max(len(seq) for seq in knot_times_batch)

# Pad the contact sequences to the maximum sequence length
padded_contact_sequences = pad_sequences(
    one_hot_contact_sequences, maxlen=max_sequence_length, dtype="float32", padding="post"
)

# Convert the padded contact sequences to a PyTorch tensor
padded_contact_sequences = torch.tensor(padded_contact_sequences)

q0 = torch.tensor(q0_batch)  # shape: (batch_size, state_size)
qf = torch.tensor(qf_batch)  # shape: (batch_size, state_size)

# Add an extra dimension to the 2D tensors
task_descriptions = task_descriptions.unsqueeze(1)  # shape: (batch_size, 1, embedding_size)
q0 = q0.unsqueeze(1)  # shape: (batch_size, 1, state_size)
qf = qf.unsqueeze(1)  # shape: (batch_size, 1, state_size)

# Expand the 2D tensors to match the size of the 3D tensors
task_descriptions = task_descriptions.expand(-1, knot_times.size(1), -1)
q0 = q0.expand(-1, knot_times.size(1), -1)
qf = qf.expand(-1, knot_times.size(1), -1)

# Add an extra dimension to the 2D tensors
knot_times = knot_times.unsqueeze(-1)  # shape: (batch_size, max_sequence_length, 1)

# Ensure padded_contact_sequences is a 3D tensor
if len(padded_contact_sequences.shape) == 2:
    padded_contact_sequences = padded_contact_sequences.unsqueeze(-1)  # shape: (batch_size, max_sequence_length, 1)

# Now, all tensors should have 3 dimensions and can be concatenated
input_data = torch.cat([task_descriptions, knot_times, padded_contact_sequences, q0, qf], dim=-1) # shape: (batch_size, max_sequence_length, embedding_size + num_features + 2 * state_size)
# Now, input_data is your final input to the neural network

print(input_data)