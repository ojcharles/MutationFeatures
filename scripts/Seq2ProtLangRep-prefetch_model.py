##### setup
model_file_dir="/tools"
seq_path = model_file_dir+"/protT5/example_seqs.fasta"
per_residue = True 
per_residue_path = model_file_dir+"/protT5/output/per_residue_embeddings.h5" # where to store the embeddings
per_protein = True
per_protein_path = model_file_dir+"/protT5/output/per_protein_embeddings.h5" # where to store the embeddings
sec_struct = True
sec_struct_path = model_file_dir+"/protT5/output/ss3_preds.fasta" # file for storing predictions
assert per_protein is True or per_residue is True or sec_struct is True, print(
    "Minimally, you need to active per_residue, per_protein or sec_struct. (or any combination)")


#@title Import dependencies and check whether GPU is available. { display-mode: "form" }
from transformers import T5EncoderModel, T5Tokenizer
import torch
import h5py
import time
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print("Using {}".format(device))



def get_T5_model():
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    model = model.to(device) # move model to GPU
    model = model.eval() # set model to evaluation model
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)
    return model, tokenizer

# pre-
get_T5_model()