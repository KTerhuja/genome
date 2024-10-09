import pysam
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import argparse
from tqdm import tqdm


def identify_species(sequence):
    # Perform BLAST search
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

    # Parse the BLAST result
    blast_record = NCBIXML.read(result_handle)

    # Extract the top hit
    if blast_record.alignments:
        top_hit = blast_record.alignments[0]
        title = top_hit.title

        # Extract genus and species from the title

        parts = title.split()
        if len(parts) >= 2:
            genus = parts[1]
            species = parts[2]
            return (f"{genus} {species}", sequence)
        else:
            return "Unable to extract genus and species from top hit"
    else:
        return "No significant matches found"
    


def extract_data(input_path, output_path):

    sequences = {'name':[], 'sequence':[]}

    with pysam.FastxFile(input_path) as fastq:
        for entry in fastq:

            sequences['name'].append(entry.name)
            sequences['sequence'].append(entry.sequence)

    
    
    sub_array = sequences['sequence']


    with ThreadPoolExecutor() as executor:
    # Submit tasks individually
        futures = [executor.submit(identify_species, sub_array) for elem in sub_array]
        
        # Use tqdm to track progress as tasks complete
        results = []
        for future in tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())

    df = pd.DataFrame(results)

    df.to_csv(output_path,index=False)

    


def parse_arguments():
    # create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-path", type=str, required=True)
    parser.add_argument("--output-path", type=str, required=True)
    return parser.parse_args()


if __name__ == "__main__": 

    args = parse_arguments()

    extract_data(args.input_path, args.output_path)



